# =============================================================================
# local_boxplot_stemness_markers.R
# Boxplot facettati (scala libera) per marker stemness + neuroendocrini
# Confronto tumore primario vs linee 2D
#
# Eseguire in locale (Mac):
#   Rscript local_boxplot_stemness_markers.R
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})
select <- dplyr::select

# -----------------------------------------------------------------------------
# 1. DATI
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento dati ---")

tpm_mat <- read_csv(file.path(RESULTS_BULK, "tpm_matrix.csv"),
                    show_col_types = FALSE) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

# Usa il metadata ordinato (sample_id allineati alle colonne TPM)
meta <- read_csv(file.path(RESULTS_BULK, "sample_metadata_ordered.csv"),
                 show_col_types = FALSE)
meta$condition <- factor(meta$condition, levels = c("tumor", "2D"))

degs_full <- read_csv(file.path(RESULTS_DEGS, "degs_full_table.csv"),
                      show_col_types = FALSE)

# Strip ENSEMBL version
rownames(tpm_mat) <- sub("\\..*", "", rownames(tpm_mat))

# Mapping ENSEMBL -> symbol
sym2ens <- setNames(degs_full$gene_id, degs_full$gene_symbol)

# -----------------------------------------------------------------------------
# 2. GENI TARGET
# -----------------------------------------------------------------------------
target_genes <- c(
  "SOX2", "PROM1", "ALDH1A1", "LGR5", "EPCAM", "CD44", "KLF4",
  "BMI1", "EZH2", "MYCN", "HMGA2", "IGF2BP2", "IGF2BP3",
  "CHGA", "CHGB", "INSM1", "SYP", "ASCL1", "NEUROD1",
  "DLL3", "CD276"   # DLL3 = DDL3 corretto; CD276 = B7-H3
)

# Trova geni disponibili
available <- target_genes[target_genes %in% names(sym2ens) &
                            sym2ens[target_genes] %in% rownames(tpm_mat)]
missing   <- target_genes[!target_genes %in% available]

cat(sprintf("Geni trovati: %d/%d\n", length(available), length(target_genes)))
if (length(missing) > 0) cat("Mancanti:", paste(missing, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 3. COSTRUZIONE DATAFRAME LONG
# -----------------------------------------------------------------------------
message("\n--- 2. Preparazione dati ---")

long_df <- do.call(rbind, lapply(available, function(gene) {
  ens_id <- sym2ens[gene]
  data.frame(
    gene      = gene,
    sample_id = colnames(tpm_mat),
    tpm       = as.numeric(tpm_mat[ens_id, ]),
    stringsAsFactors = FALSE
  )
}))

long_df <- left_join(long_df, select(meta, sample_id, condition, patient),
                     by = "sample_id")
long_df$log2TPM <- log2(long_df$tpm + 1)
long_df$gene    <- factor(long_df$gene, levels = available)

# -----------------------------------------------------------------------------
# 4. BOXPLOT FACETTATO
# -----------------------------------------------------------------------------
message("\n--- 3. Creazione boxplot ---")

OUTPUT_DIR <- file.path(RESULTS_DIR, "boxplot_markers")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(long_df,
            aes(x = condition, y = log2TPM,
                color = condition, fill = condition)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.8) +
  scale_color_manual(values = COLORS_CONDITION) +
  scale_fill_manual(values  = COLORS_CONDITION) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(
    title = "Stemness & neuroendocrine markers",
    x     = NULL,
    y     = expression(log[2](TPM + 1))
  ) +
  THEME_PGL +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold.italic", size = 10))

ggsave(file.path(OUTPUT_DIR, "boxplot_stemness_neuroendocrine_markers.pdf"),
       p, width = 16, height = ceiling(length(available) / 4) * 4)

# -----------------------------------------------------------------------------
# 5. TABELLA RIASSUNTIVA
# -----------------------------------------------------------------------------
message("\n--- 4. Tabella riassuntiva ---")

summary_df <- long_df %>%
  group_by(gene, condition) %>%
  summarise(mean_log2TPM = round(mean(log2TPM), 3), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_log2TPM,
              names_prefix = "mean_log2TPM_") %>%
  left_join(
    degs_full %>%
      filter(gene_symbol %in% available) %>%
      select(gene = gene_symbol, log2FC = log2FoldChange, padj) %>%
      distinct(gene, .keep_all = TRUE),
    by = "gene"
  )

print(as.data.frame(summary_df))
write_csv(summary_df, file.path(OUTPUT_DIR, "stemness_neuroendocrine_markers_summary.csv"))

message("\n=== local_boxplot_stemness_markers.R completato ===")
message("Output in: ", OUTPUT_DIR)
