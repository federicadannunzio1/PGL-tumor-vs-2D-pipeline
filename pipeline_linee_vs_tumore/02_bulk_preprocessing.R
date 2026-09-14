# =============================================================================
# 02_bulk_preprocessing.R
# Caricamento e preprocessing bulk RNA-seq (tumore + linee primarie 2D)
# Input:  file Salmon quant.genes.sf (16 campioni)
#         data/sample_metadata.csv
# Output: counts_matrix.RDS, tpm_matrix.csv, QC plots
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# 1. CARICAMENTO METADATA
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento metadata ---")

meta <- read_csv(SAMPLE_METADATA, show_col_types = FALSE)
meta$condition <- factor(meta$condition, levels = c("tumor", "2D"))
meta$patient   <- factor(meta$patient)

cat(sprintf("Campioni totali: %d\n", nrow(meta)))
cat(sprintf("  Tumore: %d | Linee 2D: %d\n",
    sum(meta$condition == "tumor"), sum(meta$condition == "2D")))
cat(sprintf("  Coppie paired: %d\n", sum(meta$paired)))

# -----------------------------------------------------------------------------
# 2. COSTRUZIONE PERCORSI FILE SALMON
# -----------------------------------------------------------------------------
message("\n--- 2. Ricerca file Salmon ---")

salmon_files <- file.path(BULK_SALMON_DIR, meta$file_gene)
names(salmon_files) <- meta$sample_id

missing_files <- salmon_files[!file.exists(salmon_files)]
if (length(missing_files) > 0) {
  stop("File Salmon mancanti:\n", paste(names(missing_files), collapse = "\n"))
}
message(sprintf("Trovati %d file Salmon.", length(salmon_files)))

# -----------------------------------------------------------------------------
# 3. IMPORT CON TXIMPORT (gene-level)
# ignoreTxVersion = TRUE: rimuove il suffisso di versione dagli ENSEMBL ID
# (es. ENSG00000000003.15 -> ENSG00000000003)
# -----------------------------------------------------------------------------
message("\n--- 3. Import con tximport ---")

txi <- tximport(
  files          = salmon_files,
  type           = "salmon",
  txIn           = FALSE,      # file gia' gene-level (quant.genes.sf)
  txOut          = FALSE,
  ignoreTxVersion = TRUE
)

cat(sprintf("Geni importati: %d\n", nrow(txi$counts)))
cat(sprintf("Campioni: %d\n", ncol(txi$counts)))

# -----------------------------------------------------------------------------
# 4. QC PRE-FILTRO
# -----------------------------------------------------------------------------
message("\n--- 4. QC ---")

# Totale read per campione
total_counts <- colSums(txi$counts)
qc_df <- data.frame(
  sample_id = names(total_counts),
  total_counts = total_counts,
  condition = meta$condition[match(names(total_counts), meta$sample_id)]
)

p_totalcounts <- ggplot(qc_df, aes(x = reorder(sample_id, total_counts),
                                    y = total_counts / 1e6,
                                    fill = condition)) +
  geom_col() +
  scale_fill_manual(values = COLORS_CONDITION) +
  coord_flip() +
  labs(title = "Total reads per sample",
       x = NULL, y = "Millions of reads", fill = "Condition") +
  THEME_PGL

ggsave(file.path(RESULTS_BULK, "qc_total_counts.pdf"),
       p_totalcounts, width = 8, height = 6)

# Geni espressi per campione (TPM > 1)
n_expressed <- colSums(txi$abundance > 1)
qc_df$n_expressed <- n_expressed[match(qc_df$sample_id, names(n_expressed))]

p_expressed <- ggplot(qc_df, aes(x = reorder(sample_id, n_expressed),
                                   y = n_expressed,
                                   fill = condition)) +
  geom_col() +
  scale_fill_manual(values = COLORS_CONDITION) +
  coord_flip() +
  labs(title = "Expressed genes per sample (TPM > 1)",
       x = NULL, y = "No. of genes", fill = "Condition") +
  THEME_PGL

ggsave(file.path(RESULTS_BULK, "qc_genes_expressed.pdf"),
       p_expressed, width = 8, height = 6)

# -----------------------------------------------------------------------------
# 5. PCA ESPLORATIVA (su log2 TPM)
# -----------------------------------------------------------------------------
message("\n--- 5. PCA esplorativa ---")

# Filtro geni con TPM medio > 1 per la PCA
tpm_mat <- txi$abundance
tpm_filt <- tpm_mat[rowMeans(tpm_mat) > 1, ]
log_tpm  <- log2(tpm_filt + 1)

pca_res  <- prcomp(t(log_tpm), scale. = TRUE)
pca_df   <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample_id <- rownames(pca_df)
pca_df   <- left_join(pca_df, meta, by = "sample_id")

var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                              color = condition, label = patient)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = COLORS_CONDITION) +
  labs(
    title = "PCA - Bulk RNA-seq (log2 TPM)",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2]),
    color = "Condition"
  ) +
  THEME_PGL

ggsave(file.path(RESULTS_BULK, "qc_pca.pdf"), p_pca, width = 7, height = 6)

# -----------------------------------------------------------------------------
# 6. HEATMAP CORRELAZIONE
# -----------------------------------------------------------------------------
message("\n--- 6. Heatmap correlazione ---")

cor_mat <- cor(log_tpm, method = "pearson")

# Annotazione colonne
ann_df <- data.frame(
  row.names = colnames(cor_mat),
  Condition = meta$condition[match(colnames(cor_mat), meta$sample_id)],
  Patient   = meta$patient[match(colnames(cor_mat), meta$sample_id)]
)
ann_colors <- list(Condition = COLORS_CONDITION)

pdf(file.path(RESULTS_BULK, "qc_correlation_heatmap.pdf"), width = 10, height = 9)
pheatmap(
  cor_mat,
  annotation_col  = ann_df,
  annotation_colors = ann_colors,
  color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks          = seq(0.7, 1.0, length.out = 101),
  main            = "Pearson Correlation (log2 TPM)",
  fontsize        = 9,
  border_color    = NA
)
dev.off()

# -----------------------------------------------------------------------------
# 7. SALVATAGGIO OUTPUT
# -----------------------------------------------------------------------------
message("\n--- 7. Salvataggio output ---")

# Oggetto tximport completo (per DESeq2)
saveRDS(txi,  file.path(RESULTS_BULK, "tximport_object.RDS"))

# Matrice TPM (per deconvoluzione)
write_csv(
  as.data.frame(tpm_mat) %>% tibble::rownames_to_column("gene_id"),
  file.path(RESULTS_BULK, "tpm_matrix.csv")
)

# Metadata allineato con l'ordine delle colonne di txi
write_csv(
  meta[match(colnames(txi$counts), meta$sample_id), ],
  file.path(RESULTS_BULK, "sample_metadata_ordered.csv")
)

message("\n=== 02_bulk_preprocessing.R completato ===")
message("Output in: ", RESULTS_BULK)
