# =============================================================================
# 04_deconvolution.R
# Deconvoluzione bulk RNA-seq (tumore + linee 2D) con MuSiC
# usando l'scRNA-seq del tumore come reference
#
# Obiettivo: mostrare che le linee 2D NON sono composte al 100% da cellule
# tumorali e che contengono cellule mesenchimali, come nel tumore originale
#
# Input:  results/scrna_reference/scrna_expressionset_for_music.RDS
#         results/bulk_preprocessing/tpm_matrix.csv
#         data/sample_metadata.csv
# Output: proporzioni cellulari per campione, figure comparative
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(MuSiC)
  library(Biobase)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

check_inputs(
  file.path(RESULTS_SCRNA, "scrna_expressionset_for_music.RDS"),
  file.path(RESULTS_BULK,  "tpm_matrix.csv"),
  SAMPLE_METADATA
)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO DATI
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento dati ---")

scrna_eset <- readRDS(
  file.path(RESULTS_SCRNA, "scrna_expressionset_for_music.RDS")
)

# MuSiC v1.0.0 requires SingleCellExperiment (not ExpressionSet)
scrna_sce <- SingleCellExperiment(
  assays  = list(counts = Biobase::exprs(scrna_eset)),
  colData = Biobase::pData(scrna_eset)
)

cat(sprintf("Reference scRNA: %d geni x %d cellule\n",
            nrow(assay(scrna_sce)), ncol(assay(scrna_sce))))
cat("Tipi cellulari nel reference:\n")
print(table(scrna_sce$cellType))

# MuSiC richiede >= 2 soggetti nel reference per stimare la varianza.
# Se il reference ha un solo sampleID (tumore singolo), creiamo pseudo-repliche
# dividendo casualmente le cellule di ogni tipo in 3 gruppi.
n_subjects <- length(unique(scrna_sce$sampleID))
cat(sprintf("\nSoggetti unici nel reference scRNA: %d\n", n_subjects))
if (n_subjects < 2) {
  message("Reference ha 1 solo soggetto: creazione di 3 pseudo-repliche per MuSiC...")
  set.seed(SEED)
  ct_vec  <- scrna_sce$cellType
  ps_vec  <- character(ncol(scrna_sce))
  for (ct in unique(ct_vec)) {
    idx <- which(ct_vec == ct)
    ps_vec[idx] <- sample(paste0("subj_", 1:3), length(idx), replace = TRUE)
  }
  scrna_sce$sampleID <- ps_vec
  cat("Distribuzione pseudo-soggetti:\n")
  print(table(scrna_sce$sampleID))
}

tpm_mat <- read_csv(file.path(RESULTS_BULK, "tpm_matrix.csv"),
                    show_col_types = FALSE) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

meta <- read_csv(SAMPLE_METADATA, show_col_types = FALSE)
meta$condition <- factor(meta$condition, levels = c("tumor", "2D"))

cat(sprintf("\nBulk RNA-seq: %d geni x %d campioni\n",
            nrow(tpm_mat), ncol(tpm_mat)))

# Allinea nomi campioni tra metadata e colonne TPM
tpm_mat <- tpm_mat[, meta$sample_id]

# -----------------------------------------------------------------------------
# 1b. VERIFICA OVERLAP GENI bulk vs scRNA (e conversione se necessario)
# Caso comune: bulk usa ENSEMBL IDs, scRNA usa simboli genici
# -----------------------------------------------------------------------------
message("\n--- 1b. Verifica overlap geni ---")

common_genes_pre <- intersect(rownames(tpm_mat), rownames(scrna_sce))
cat(sprintf("Geni in comune bulk vs scRNA (prima di eventuale conversione): %d\n",
            length(common_genes_pre)))
cat("Esempi bulk gene IDs: ",  paste(head(rownames(tpm_mat),  4), collapse = ", "), "\n")
cat("Esempi scRNA gene IDs:", paste(head(rownames(scrna_sce), 4), collapse = ", "), "\n")

if (length(common_genes_pre) < 100) {
  bulk_is_ensembl  <- grepl("^ENSG", rownames(tpm_mat)[1])
  scrna_is_ensembl <- grepl("^ENSG", rownames(scrna_sce)[1])

  if (bulk_is_ensembl && !scrna_is_ensembl) {
    message("Bulk usa ENSEMBL, scRNA usa simboli -> conversione ENSEMBL->simbolo sul bulk...")
    # Rimuovi numeri di versione (es. ENSG00000223972.5 -> ENSG00000223972)
    rownames(tpm_mat) <- sub("\\..*", "", rownames(tpm_mat))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    emap <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = rownames(tpm_mat),
      columns = "SYMBOL",
      keytype = "ENSEMBL"
    )
    emap <- emap[!is.na(emap$SYMBOL) & !duplicated(emap$ENSEMBL), ]
    idx  <- match(rownames(tpm_mat), emap$ENSEMBL)
    keep <- !is.na(idx)
    tpm_mat2 <- tpm_mat[keep, , drop = FALSE]
    rownames(tpm_mat2) <- emap$SYMBOL[idx[keep]]
    tpm_mat  <- tpm_mat2[!duplicated(rownames(tpm_mat2)), , drop = FALSE]
    cat(sprintf("Bulk dopo conversione: %d geni\n", nrow(tpm_mat)))
    common_genes_post <- intersect(rownames(tpm_mat), rownames(scrna_sce))
    cat(sprintf("Geni in comune dopo conversione: %d\n", length(common_genes_post)))
    if (length(common_genes_post) < 100)
      stop("Ancora troppo pochi geni in comune dopo conversione ENSEMBL->simbolo.")
  } else if (!bulk_is_ensembl && scrna_is_ensembl) {
    message("scRNA usa ENSEMBL, bulk usa simboli -> conversione ENSEMBL->simbolo su scRNA...")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    emap <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = rownames(scrna_sce),
      columns = "SYMBOL",
      keytype = "ENSEMBL"
    )
    emap <- emap[!is.na(emap$SYMBOL) & !duplicated(emap$ENSEMBL), ]
    idx  <- match(rownames(scrna_sce), emap$ENSEMBL)
    keep <- !is.na(idx)
    sce_counts2 <- assay(scrna_sce, "counts")[keep, , drop = FALSE]
    rownames(sce_counts2) <- emap$SYMBOL[idx[keep]]
    sce_counts2 <- sce_counts2[!duplicated(rownames(sce_counts2)), , drop = FALSE]
    scrna_sce <- SingleCellExperiment(
      assays  = list(counts = sce_counts2),
      colData = colData(scrna_sce)
    )
    cat(sprintf("scRNA dopo conversione: %d geni\n", nrow(scrna_sce)))
    common_genes_post <- intersect(rownames(tpm_mat), rownames(scrna_sce))
    cat(sprintf("Geni in comune dopo conversione: %d\n", length(common_genes_post)))
    if (length(common_genes_post) < 100)
      stop("Ancora troppo pochi geni in comune dopo conversione ENSEMBL->simbolo.")
  } else {
    stop(sprintf(
      "Overlap insufficiente (%d geni) e formato identico. Controllare i dati.",
      length(common_genes_pre)
    ))
  }
} else {
  cat("Overlap sufficiente - nessuna conversione necessaria.\n")
}

# -----------------------------------------------------------------------------
# 2. COSTRUZIONE ExpressionSet PER IL BULK
# MuSiC richiede che il bulk sia un ExpressionSet
# -----------------------------------------------------------------------------
message("\n--- 2. Costruzione ExpressionSet bulk ---")

bulk_eset <- ExpressionSet(
  assayData = tpm_mat,
  phenoData = new("AnnotatedDataFrame",
                  data = data.frame(
                    row.names  = colnames(tpm_mat),
                    sample_id  = meta$sample_id,
                    condition  = meta$condition,
                    patient    = meta$patient
                  ))
)

# -----------------------------------------------------------------------------
# 3. DECONVOLUZIONE CON MuSiC
# music_prop(): stima le proporzioni di tipi cellulari per ogni campione bulk
# usando come reference le medie pesate dal reference scRNA-seq multi-soggetto
# -----------------------------------------------------------------------------
message("\n--- 3. Deconvoluzione MuSiC ---")
message("Questo passo puo' richiedere alcuni minuti...")

set.seed(SEED)

music_results <- music_prop(
  bulk.mtx    = exprs(bulk_eset),
  sc.sce      = scrna_sce,
  clusters    = "cellType",
  samples     = "sampleID",
  select.ct   = NULL,              # usa tutti i tipi cellulari
  verbose     = FALSE
)

# Estrai la matrice delle proporzioni (campioni x tipi cellulari)
prop_mat <- music_results$Est.prop.weighted

cat("\nProporzioni stimate (prime righe):\n")
print(round(head(prop_mat, 5), 3))

# Controlla: le proporzioni sommano a 1?
cat("\nSomma proporzioni per campione (dovrebbero essere ~1):\n")
print(round(rowSums(prop_mat), 3))

# -----------------------------------------------------------------------------
# 4. FORMATTAZIONE RISULTATI
# -----------------------------------------------------------------------------
message("\n--- 4. Formattazione risultati ---")

prop_df <- as.data.frame(prop_mat) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(select(meta, sample_id, condition, patient, paired),
            by = "sample_id") %>%
  pivot_longer(
    cols      = -c(sample_id, condition, patient, paired),
    names_to  = "cell_type",
    values_to = "proportion"
  )

# Flag mesenchimali
mesenchymal_pattern <- "(?i)(mesench|fibroblast|stromal|sustentacular|MSC)"
prop_df$is_mesenchymal <- grepl(mesenchymal_pattern,
                                 prop_df$cell_type, perl = TRUE)

write_csv(prop_df, file.path(RESULTS_DECONV, "music_proportions_long.csv"))

# Formato wide per lettura facile
prop_wide <- as.data.frame(prop_mat) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(select(meta, sample_id, condition, patient), by = "sample_id") %>%
  relocate(condition, patient, .after = sample_id)

write_csv(prop_wide, file.path(RESULTS_DECONV, "music_proportions_wide.csv"))

cat("\nRiepilogo proporzione mesenchimale:\n")
mes_summary <- prop_df %>%
  filter(is_mesenchymal) %>%
  group_by(sample_id, condition, patient) %>%
  summarise(prop_mesenchymal = sum(proportion), .groups = "drop") %>%
  arrange(condition, desc(prop_mesenchymal))
print(mes_summary)
write_csv(mes_summary, file.path(RESULTS_DECONV, "mesenchymal_proportions_bulk.csv"))

# -----------------------------------------------------------------------------
# 5. TEST STATISTICO: confronto proporzioni tumore vs 2D
# Wilcoxon paired test su campioni con coppie complete (base R)
# -----------------------------------------------------------------------------
message("\n--- 5. Test statistici (Wilcoxon paired) ---")

# Solo coppie paired
prop_paired <- prop_df %>% filter(paired == TRUE)

cell_types_all <- unique(prop_paired$cell_type)

stat_results <- do.call(rbind, lapply(cell_types_all, function(ct) {
  df_ct <- prop_paired %>% filter(cell_type == ct) %>%
    arrange(patient, condition)
  tumor_vals <- df_ct %>% filter(condition == "tumor") %>% pull(proportion)
  line_vals  <- df_ct %>% filter(condition == "2D")    %>% pull(proportion)

  # Wilcoxon paired (gestisce anche coppie incomplete)
  res <- tryCatch(
    wilcox.test(tumor_vals, line_vals, paired = TRUE, exact = FALSE),
    error = function(e) list(statistic = NA, p.value = NA)
  )
  data.frame(
    cell_type = ct,
    statistic = as.numeric(res$statistic),
    p         = res$p.value,
    stringsAsFactors = FALSE
  )
}))

stat_results$p.adj <- p.adjust(stat_results$p, method = "BH")
stat_results$significance <- cut(
  stat_results$p.adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right  = TRUE
)
stat_results <- stat_results[order(stat_results$p.adj), ]

cat("\nTest Wilcoxon paired per tipo cellulare:\n")
print(stat_results)
write_csv(stat_results,
          file.path(RESULTS_DECONV, "wilcoxon_celltype_tumor_vs_2D.csv"))

# Tipi cellulari significativamente diversi
sig_types <- stat_results$cell_type[!is.na(stat_results$p.adj) &
                                      stat_results$p.adj < 0.05]

cat(sprintf("\nTipi cellulari significativamente diversi (padj < 0.05): %d\n",
            length(sig_types)))
if (length(sig_types) > 0) print(sig_types)

# -----------------------------------------------------------------------------
# 6. FIGURE
# -----------------------------------------------------------------------------
message("\n--- 6. Figure ---")

# Palette colori per tipi cellulari
n_types     <- length(unique(prop_df$cell_type))
type_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_types),
  unique(prop_df$cell_type)
)

## 6a. Stacked barplot - composizione per campione
# Campioni ordinati per paziente, tumore prima poi 2D
sample_order <- meta %>%
  arrange(patient, condition) %>%
  pull(sample_id)

prop_df$sample_id <- factor(prop_df$sample_id, levels = sample_order)
prop_df$condition <- factor(prop_df$condition, levels = c("tumor", "2D"))

p_stack <- ggplot(prop_df,
                  aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = type_colors) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Estimated cell type composition (MuSiC)",
    subtitle = "Bulk RNA-seq deconvolution using tumor scRNA-seq as reference",
    x = NULL, y = "Estimated proportion", fill = "Cell type"
  ) +
  THEME_PGL +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 8),
    legend.text   = element_text(size = 8),
    strip.text    = element_text(face = "bold", size = 11)
  )

ggsave(file.path(RESULTS_DECONV_FIG, "stacked_barplot_proportions.pdf"),
       p_stack, width = 14, height = 7)

## 6b. Boxplot paired: proporzione mesenchimale tumore vs 2D
mes_box_df <- prop_df %>%
  filter(is_mesenchymal) %>%
  group_by(sample_id, condition, patient) %>%
  summarise(prop_mes = sum(proportion), .groups = "drop")

# Calcola p-value manuale per mesenchimale (Wilcoxon paired)
mes_tumor <- mes_box_df %>% filter(condition == "tumor") %>%
  arrange(patient) %>% pull(prop_mes)
mes_2d    <- mes_box_df %>% filter(condition == "2D") %>%
  arrange(patient) %>% pull(prop_mes)
mes_pval  <- tryCatch(
  wilcox.test(mes_tumor, mes_2d, paired = TRUE, exact = FALSE)$p.value,
  error = function(e) NA
)
mes_plab  <- ifelse(is.na(mes_pval), "p = NA",
                    ifelse(mes_pval < 0.001, "p < 0.001",
                           sprintf("p = %.3f", mes_pval)))
mes_ymax  <- max(mes_box_df$prop_mes * 100, na.rm = TRUE)

p_mes_box <- ggplot(mes_box_df,
                    aes(x = condition, y = prop_mes * 100,
                        color = condition, fill = condition)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(aes(group = patient), color = "grey40", alpha = 0.6, linewidth = 0.7) +
  scale_color_manual(values = COLORS_CONDITION) +
  scale_fill_manual(values  = COLORS_CONDITION) +
  annotate("segment",
           x = 1, xend = 2,
           y = mes_ymax * 1.08, yend = mes_ymax * 1.08,
           color = "black") +
  annotate("text",
           x = 1.5, y = mes_ymax * 1.12,
           label = mes_plab, size = 3.5) +
  labs(
    title    = "Mesenchymal cell proportion",
    subtitle = "Primary tumor vs 2D cell lines (paired Wilcoxon test)",
    x = "Condition", y = "% mesenchymal cells",
    color = NULL, fill = NULL
  ) +
  THEME_PGL + theme(legend.position = "none")

ggsave(file.path(RESULTS_DECONV_FIG, "boxplot_mesenchymal_tumor_vs_2D.pdf"),
       p_mes_box, width = 6, height = 7)

## 6c. Boxplot paired per ogni tipo cellulare significativo
if (length(sig_types) > 0) {
  prop_sig <- prop_df %>%
    filter(cell_type %in% sig_types, paired == TRUE)

  p_sig_types <- ggplot(prop_sig,
                        aes(x = condition, y = proportion * 100,
                            color = condition, fill = condition)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    geom_point(size = 2) +
    geom_line(aes(group = patient), color = "grey40", alpha = 0.5,
              linewidth = 0.5) +
    scale_color_manual(values = COLORS_CONDITION) +
    scale_fill_manual(values  = COLORS_CONDITION) +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(
      title = "Significantly different cell types (padj < 0.05)",
      x = "Condition", y = "% estimated proportion",
      color = NULL, fill = NULL
    ) +
    THEME_PGL + theme(legend.position = "none")

  ggsave(file.path(RESULTS_DECONV_FIG, "boxplot_significant_celltypes.pdf"),
         p_sig_types,
         width  = min(4 * ceiling(sqrt(length(sig_types))), 16),
         height = min(4 * ceiling(length(sig_types) /
                        ceiling(sqrt(length(sig_types)))), 16))
}

## 6d. Heatmap proporzioni (campioni x tipi cellulari)
prop_heatmap <- prop_mat

# Ordina campioni come stacked barplot
prop_heatmap <- prop_heatmap[sample_order, ]

ann_row <- data.frame(
  row.names = rownames(prop_heatmap),
  Condition = meta$condition[match(rownames(prop_heatmap), meta$sample_id)],
  Patient   = meta$patient[match(rownames(prop_heatmap), meta$sample_id)]
)

pdf(file.path(RESULTS_DECONV_FIG, "heatmap_proportions.pdf"),
    width = 10, height = 8)
pheatmap::pheatmap(
  prop_heatmap,
  annotation_row  = ann_row,
  annotation_colors = list(Condition = COLORS_CONDITION),
  color           = colorRampPalette(c("white", "#E64B35"))(100),
  cluster_rows    = FALSE,
  cluster_cols    = TRUE,
  main            = "Estimated cell type proportions (MuSiC)",
  fontsize         = 9,
  border_color    = NA
)
dev.off()

# -----------------------------------------------------------------------------
# 7. SALVATAGGIO OGGETTO MuSiC COMPLETO
# -----------------------------------------------------------------------------
message("\n--- 7. Salvataggio ---")

saveRDS(music_results,
        file.path(RESULTS_DECONV, "music_results_full.RDS"))

message("\n=== 04_deconvolution.R completato ===")
message("Output in: ", RESULTS_DECONV)
