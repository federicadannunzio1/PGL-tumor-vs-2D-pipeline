# =============================================================================
# 03_degs.R
# DESeq2: analisi differenziale tumore vs linee primarie 2D (bulk RNA-seq)
# Input:  results/bulk_preprocessing/tximport_object.RDS
#         data/sample_metadata.csv
# Output: tabella DEG completa, filtrata, arricchimento GO/KEGG, figure
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# clusterProfiler opzionale: GSEA viene saltato se non installato
HAS_CLUSTERPROFILER <- requireNamespace("clusterProfiler", quietly = TRUE)
if (HAS_CLUSTERPROFILER) {
  suppressPackageStartupMessages(library(clusterProfiler))
  message("clusterProfiler disponibile: GSEA verra' eseguita.")
} else {
  message("clusterProfiler non disponibile: la sezione GSEA verra' saltata.")
}

check_inputs(
  file.path(RESULTS_BULK, "tximport_object.RDS"),
  SAMPLE_METADATA
)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO DATI
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento dati ---")

txi  <- readRDS(file.path(RESULTS_BULK, "tximport_object.RDS"))
meta <- read_csv(SAMPLE_METADATA, show_col_types = FALSE)
meta$condition <- factor(meta$condition, levels = c("tumor", "2D"))
meta$patient   <- factor(meta$patient)

# Allinea metadata con ordine colonne txi
meta <- meta[match(colnames(txi$counts), meta$sample_id), ]
stopifnot(all(colnames(txi$counts) == meta$sample_id))

# -----------------------------------------------------------------------------
# 2. COSTRUZIONE DESeqDataSet
# Strategia: disegno paired usando solo le 7 coppie complete
# ~ patient + condition  (patient come blocking factor)
# -----------------------------------------------------------------------------
message("\n--- 2. Costruzione DESeqDataSet ---")

# Subset: solo campioni paired
meta_paired <- meta[meta$paired == TRUE, ]
txi_counts_paired <- txi$counts[, meta_paired$sample_id]
txi_paired <- list(
  counts    = txi$counts[, meta_paired$sample_id],
  abundance = txi$abundance[, meta_paired$sample_id],
  length    = txi$length[, meta_paired$sample_id]
)
class(txi_paired) <- "list"

cat(sprintf("Campioni nel disegno paired: %d (%d pazienti x 2)\n",
            nrow(meta_paired), sum(meta_paired$condition == "tumor")))

dds <- DESeqDataSetFromTximport(
  txi     = txi_paired,
  colData = meta_paired,
  design  = ~ patient + condition
)

# -----------------------------------------------------------------------------
# 3. FILTRO GENI A BASSA ESPRESSIONE
# -----------------------------------------------------------------------------
message("\n--- 3. Filtro geni a bassa espressione ---")

keep <- rowSums(counts(dds) >= DEG_MIN_COUNTS) >= DEG_MIN_SAMPLES
dds  <- dds[keep, ]
cat(sprintf("Geni prima del filtro: %d\n", sum(!keep) + sum(keep)))
cat(sprintf("Geni dopo filtro (>= %d counts in >= %d campioni): %d\n",
            DEG_MIN_COUNTS, DEG_MIN_SAMPLES, sum(keep)))

# -----------------------------------------------------------------------------
# 4. NORMALIZZAZIONE E DESeq2
# -----------------------------------------------------------------------------
message("\n--- 4. DESeq2 ---")

dds <- DESeq(dds, parallel = (N_CORES > 1))

# Dispersion plot
pdf(file.path(RESULTS_DEGS_FIG, "dispersion_plot.pdf"), width = 7, height = 5)
plotDispEsts(dds, main = "Dispersion estimates - DESeq2")
dev.off()

# -----------------------------------------------------------------------------
# 5. ESTRAZIONE RISULTATI CON LFC SHRINKAGE
# Confronto: 2D vs tumor (reference = tumor)
# Positivo log2FC = piu' espresso nelle linee 2D
# Negativo log2FC = piu' espresso nel tumore
# -----------------------------------------------------------------------------
message("\n--- 5. Risultati con lfcShrink (apeglm) ---")

res_lfc <- lfcShrink(
  dds,
  coef     = "condition_2D_vs_tumor",
  type     = "apeglm",
  parallel = (N_CORES > 1)
)

# Anche risultati non shrinkati (per il volcano plot)
res_raw <- results(dds, contrast = c("condition", "2D", "tumor"),
                   independentFiltering = TRUE)

cat(sprintf("\nRiepilogo risultati (soglie: padj < %.2f, |LFC| > %.1f):\n",
            DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD))
summary_df <- data.frame(
  Totale_geni_testati = nrow(res_lfc),
  DEG_totali = sum(res_lfc$padj < DEG_PADJ_THRESHOLD &
                   abs(res_lfc$log2FoldChange) > DEG_LFC_THRESHOLD,
                   na.rm = TRUE),
  Up_in_2D   = sum(res_lfc$padj < DEG_PADJ_THRESHOLD &
                   res_lfc$log2FoldChange > DEG_LFC_THRESHOLD,
                   na.rm = TRUE),
  Down_in_2D = sum(res_lfc$padj < DEG_PADJ_THRESHOLD &
                   res_lfc$log2FoldChange < -DEG_LFC_THRESHOLD,
                   na.rm = TRUE)
)
print(summary_df)
write_csv(summary_df, file.path(RESULTS_DEGS, "degs_summary.csv"))

# -----------------------------------------------------------------------------
# 6. ANNOTAZIONE GENI
# ENSEMBL ID (senza versione) -> gene symbol + nome
# -----------------------------------------------------------------------------
message("\n--- 6. Annotazione gene symbols ---")

res_df <- as.data.frame(res_lfc) %>%
  tibble::rownames_to_column("gene_id")

res_df$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$gene_id,
  column  = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_df$gene_name <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$gene_id,
  column  = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_df$entrez_id <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$gene_id,
  column  = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Ordina per padj
res_df <- res_df %>% arrange(padj)

# Tabella completa
write_csv(res_df, file.path(RESULTS_DEGS, "degs_full_table.csv"))

# Tabella filtrata: solo DEG significativi
degs <- res_df %>%
  filter(padj < DEG_PADJ_THRESHOLD,
         abs(log2FoldChange) > DEG_LFC_THRESHOLD)

write_csv(degs, file.path(RESULTS_DEGS, "degs_filtered.csv"))
cat(sprintf("DEG salvati: %d\n", nrow(degs)))

# -----------------------------------------------------------------------------
# 7. CHECK MARKER MESENCHIMALI E NEUROENDOCRINI NEI DEG
# -----------------------------------------------------------------------------
message("\n--- 7. Check marker mesenchimali e neuroendocrini ---")

mes_in_degs <- degs %>%
  filter(gene_symbol %in% MESENCHYMAL_MARKERS) %>%
  select(gene_id, gene_symbol, log2FoldChange, padj) %>%
  arrange(log2FoldChange)

ne_in_degs <- degs %>%
  filter(gene_symbol %in% NEUROENDOCRINE_MARKERS) %>%
  select(gene_id, gene_symbol, log2FoldChange, padj) %>%
  arrange(log2FoldChange)

cat("\nMarker MESENCHIMALI nei DEG:\n")
print(mes_in_degs)
cat("\nMarker NEUROENDOCRINI nei DEG:\n")
print(ne_in_degs)

write_csv(mes_in_degs, file.path(RESULTS_DEGS, "mesenchymal_markers_in_degs.csv"))
write_csv(ne_in_degs,  file.path(RESULTS_DEGS, "neuroendocrine_markers_in_degs.csv"))

# Tabella anche su tutti i geni (non solo DEG)
markers_all <- res_df %>%
  filter(gene_symbol %in% c(MESENCHYMAL_MARKERS, NEUROENDOCRINE_MARKERS)) %>%
  mutate(marker_type = case_when(
    gene_symbol %in% MESENCHYMAL_MARKERS    ~ "Mesenchymal",
    gene_symbol %in% NEUROENDOCRINE_MARKERS ~ "Neuroendocrine"
  )) %>%
  arrange(marker_type, log2FoldChange)

write_csv(markers_all, file.path(RESULTS_DEGS, "known_markers_all_genes.csv"))

# -----------------------------------------------------------------------------
# 8. FIGURE
# -----------------------------------------------------------------------------
message("\n--- 8. Figure ---")

## 8a. MA plot
pdf(file.path(RESULTS_DEGS_FIG, "ma_plot.pdf"), width = 7, height = 5)
plotMA(res_raw, alpha = DEG_PADJ_THRESHOLD,
       main = "MA plot - 2D vs tumor", ylim = c(-8, 8))
abline(h = c(-DEG_LFC_THRESHOLD, DEG_LFC_THRESHOLD), col = "blue", lty = 2)
dev.off()

## 8b. Volcano plot
volcano_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < DEG_PADJ_THRESHOLD & log2FoldChange >  DEG_LFC_THRESHOLD ~ "Up in 2D",
      padj < DEG_PADJ_THRESHOLD & log2FoldChange < -DEG_LFC_THRESHOLD ~ "Down in 2D (up in tumor)",
      TRUE ~ "NS"
    ),
    label = case_when(
      gene_symbol %in% c(MESENCHYMAL_MARKERS, NEUROENDOCRINE_MARKERS) ~ gene_symbol,
      TRUE ~ NA_character_
    )
  )

volcano_colors <- c(
  "Up in 2D"              = "#4DBBD5",
  "Down in 2D (up in tumor)" = "#E64B35",
  "NS"                    = "grey70"
)

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = -log10(padj),
                        color = sig, label = label)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_point(data = filter(volcano_df, !is.na(label)),
             size = 2.5, alpha = 1) +
  geom_text_repel(
    data = filter(volcano_df, !is.na(label)),
    size = 3, max.overlaps = 30, fontface = "italic",
    segment.color = "grey40"
  ) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-DEG_LFC_THRESHOLD, DEG_LFC_THRESHOLD),
             linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(DEG_PADJ_THRESHOLD),
             linetype = "dashed", color = "grey30") +
  labs(
    title = "Volcano plot: primary cell lines (2D) vs tumor",
    subtitle = "Highlighted: mesenchymal markers (blue) and neuroendocrine markers (red)",
    x = "log2 Fold Change (2D / tumor)",
    y = "-log10 (padj)",
    color = NULL
  ) +
  THEME_PGL

ggsave(file.path(RESULTS_DEGS_FIG, "volcano_plot.pdf"),
       p_volcano, width = 9, height = 7)

## 8c. Heatmap top 50 DEG (25 up + 25 down in 2D)
top_up   <- degs %>% filter(log2FoldChange > 0) %>%
            slice_min(padj, n = 25) %>% pull(gene_id)
top_down <- degs %>% filter(log2FoldChange < 0) %>%
            slice_min(padj, n = 25) %>% pull(gene_id)
top_genes <- c(top_up, top_down)

# Normalizzazione: VST
vsd <- vst(dds, blind = FALSE)
heatmap_mat <- assay(vsd)[top_genes, ]

# Rimpiazza gene_id con gene_symbol nelle righe
rownames(heatmap_mat) <- res_df$gene_symbol[
  match(rownames(heatmap_mat), res_df$gene_id)
]
rownames(heatmap_mat)[is.na(rownames(heatmap_mat))] <- top_genes[
  is.na(rownames(heatmap_mat))
]

# Scala per riga (z-score)
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

ann_col <- data.frame(
  row.names = colnames(heatmap_mat),
  Condition = meta_paired$condition,
  Patient   = meta_paired$patient
)

pdf(file.path(RESULTS_DEGS_FIG, "heatmap_top50_degs.pdf"), width = 10, height = 12)
pheatmap(
  heatmap_mat_scaled,
  annotation_col  = ann_col,
  annotation_colors = list(Condition = COLORS_CONDITION),
  color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks          = seq(-3, 3, length.out = 101),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 8,
  fontsize_col    = 9,
  main            = "Top 50 DEGs (z-score VST)\n2D vs tumor"
)
dev.off()

## 8d. Dotplot marker mesenchimali - espressione media per condizione
markers_expr <- assay(vsd)[
  res_df$gene_id[res_df$gene_symbol %in% MESENCHYMAL_MARKERS &
                 !is.na(res_df$gene_symbol)], ,
  drop = FALSE
]
rownames(markers_expr) <- res_df$gene_symbol[
  res_df$gene_id %in% rownames(markers_expr)
]

markers_df <- as.data.frame(t(markers_expr)) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(select(meta_paired, sample_id, condition, patient),
            by = "sample_id") %>%
  tidyr::pivot_longer(
    cols      = -c(sample_id, condition, patient),
    names_to  = "gene",
    values_to = "expression"
  )

p_markers <- ggplot(markers_df, aes(x = gene, y = expression,
                                     color = condition, fill = condition)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_color_manual(values = COLORS_CONDITION) +
  scale_fill_manual(values  = COLORS_CONDITION) +
  labs(
    title = "Mesenchymal marker expression (VST)",
    x = NULL, y = "Normalized expression (VST)",
    color = "Condition", fill = "Condition"
  ) +
  THEME_PGL +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

ggsave(file.path(RESULTS_DEGS_FIG, "mesenchymal_markers_expression.pdf"),
       p_markers, width = 11, height = 6)

# -----------------------------------------------------------------------------
# 9. GENE SET ENRICHMENT (GSEA) con clusterProfiler (opzionale)
# -----------------------------------------------------------------------------
message("\n--- 9. Gene Set Enrichment Analysis ---")

if (!HAS_CLUSTERPROFILER) {
  message("clusterProfiler non disponibile: sezione GSEA saltata.")
} else {

  # Vettore ordinato per stat (score GSEA)
  gene_list <- res_df %>%
    filter(!is.na(entrez_id), !is.na(stat)) %>%
    arrange(desc(stat))

  ranked_vec <- setNames(gene_list$stat, gene_list$entrez_id)

  # Rimuovi duplicati (conserva il piu' alto)
  ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]

  ## GSEA GO - Biological Process
  set.seed(SEED)
  gsea_go <- gseGO(
    geneList     = ranked_vec,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose      = FALSE
  )

  if (nrow(as.data.frame(gsea_go)) > 0) {
    write_csv(as.data.frame(gsea_go),
              file.path(RESULTS_DEGS, "gsea_go_bp_results.csv"))

    p_gsea_go <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = "GSEA - GO Biological Process") +
      THEME_PGL

    ggsave(file.path(RESULTS_DEGS_FIG, "gsea_go_dotplot.pdf"),
           p_gsea_go, width = 14, height = 8)
  } else {
    message("Nessun termine GO significativo trovato.")
  }

  ## GSEA KEGG
  set.seed(SEED)
  gsea_kegg <- gseKEGG(
    geneList      = ranked_vec,
    organism      = "hsa",
    minGSSize     = 15,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )

  if (nrow(as.data.frame(gsea_kegg)) > 0) {
    write_csv(as.data.frame(gsea_kegg),
              file.path(RESULTS_DEGS, "gsea_kegg_results.csv"))

    p_gsea_kegg <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = "GSEA - KEGG Pathways") +
      THEME_PGL

    ggsave(file.path(RESULTS_DEGS_FIG, "gsea_kegg_dotplot.pdf"),
           p_gsea_kegg, width = 14, height = 8)
  } else {
    message("Nessun pathway KEGG significativo trovato.")
  }

} # end if (HAS_CLUSTERPROFILER)

# -----------------------------------------------------------------------------
# 10. SALVATAGGIO OGGETTO DESeq2 (per script successivi)
# -----------------------------------------------------------------------------
message("\n--- 10. Salvataggio oggetti R ---")

saveRDS(dds, file.path(RESULTS_DEGS, "dds_object.RDS"))
saveRDS(res_df, file.path(RESULTS_DEGS, "deseq2_results_df.RDS"))
saveRDS(vsd, file.path(RESULTS_DEGS, "vsd_object.RDS"))

message("\n=== 03_degs.R completato ===")
message("Output in: ", RESULTS_DEGS)
