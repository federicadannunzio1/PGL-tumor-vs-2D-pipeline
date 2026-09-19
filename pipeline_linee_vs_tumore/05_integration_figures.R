# =============================================================================
# 05_integration_figures.R
# Integrazione di tutti i risultati e figure finali
#
# Obiettivi:
# 1. Mostrare cellule mesenchimali nel tumore originale (scRNA)
# 2. Mostrare cellule mesenchimali anche nelle linee 2D (deconvoluzione bulk)
# 3. Confermare che le linee 2D NON sono 100% tumore
# 4. Overlap tra DEG (linee vs tumore) e signature mesenchimale scRNA
#
# Input:  output di tutti gli script precedenti
# Output: figure integrate publication-ready
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggvenn)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
})

# Verifica che tutti gli output necessari esistano
check_inputs(
  file.path(RESULTS_SCRNA,  "celltype_proportions_per_sample.csv"),
  file.path(RESULTS_SCRNA,  "mesenchymal_proportions.csv"),
  file.path(RESULTS_SCRNA,  "markers_per_celltype.csv"),
  file.path(RESULTS_DEGS,   "degs_filtered.csv"),
  file.path(RESULTS_DEGS,   "degs_full_table.csv"),
  file.path(RESULTS_DECONV, "music_proportions_long.csv"),
  file.path(RESULTS_DECONV, "mesenchymal_proportions_bulk.csv")
)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO RISULTATI
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento risultati ---")

# scRNA - proporzioni
scrna_proportions <- read_csv(
  file.path(RESULTS_SCRNA, "celltype_proportions_per_sample.csv"),
  show_col_types = FALSE
)
scrna_mes_prop <- read_csv(
  file.path(RESULTS_SCRNA, "mesenchymal_proportions.csv"),
  show_col_types = FALSE
)

# scRNA - marker per tipo cellulare
scrna_markers <- read_csv(
  file.path(RESULTS_SCRNA, "markers_per_celltype.csv"),
  show_col_types = FALSE
)

# DEG
degs_all <- read_csv(
  file.path(RESULTS_DEGS, "degs_full_table.csv"),
  show_col_types = FALSE
)
degs_filtered <- read_csv(
  file.path(RESULTS_DEGS, "degs_filtered.csv"),
  show_col_types = FALSE
)

# Deconvoluzione bulk
bulk_prop <- read_csv(
  file.path(RESULTS_DECONV, "music_proportions_long.csv"),
  show_col_types = FALSE
)
bulk_mes_prop <- read_csv(
  file.path(RESULTS_DECONV, "mesenchymal_proportions_bulk.csv"),
  show_col_types = FALSE
)

# Marker mesenchimali da Pasquale (o calcolati)
mes_markers_in_degs <- read_csv(
  file.path(RESULTS_DEGS, "mesenchymal_markers_in_degs.csv"),
  show_col_types = FALSE
)

message("Tutti i file caricati.")

# -----------------------------------------------------------------------------
# 2. PANEL A - Cellule mesenchimali nel tumore originale (scRNA)
# -----------------------------------------------------------------------------
message("\n--- 2. Panel A: mesenchimali nel tumore (scRNA) ---")

n_types <- length(unique(scrna_proportions$cell_type))
cell_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_types),
  unique(scrna_proportions$cell_type)
)

panel_a_stack <- ggplot(
  scrna_proportions,
  aes(x = sample, y = proportion, fill = cell_type)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = cell_colors) +
  labs(
    title = "A \u2014 Tumor cell type composition",
    x = NULL, y = "Proportion", fill = "Cell type"
  ) +
  THEME_PGL +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

panel_a_mes <- ggplot(
  scrna_mes_prop,
  aes(x = reorder(sample, proportion_mesenchymal),
      y = proportion_mesenchymal * 100)
) +
  geom_col(fill = "#F39B7F", width = 0.75) +
  geom_text(
    aes(label = sprintf("%.1f%%", proportion_mesenchymal * 100)),
    hjust = -0.1, size = 3.5
  ) +
  coord_flip() +
  ylim(0, max(scrna_mes_prop$proportion_mesenchymal * 100) * 1.25) +
  labs(
    title = "A2 \u2014 Mesenchymal cells per patient",
    x = "Sample", y = "% mesenchymal cells"
  ) +
  THEME_PGL

pdf(file.path(RESULTS_FIGURES, "panelA_scrna_tumor_composition.pdf"),
    width = 16, height = 7)
print(panel_a_stack + panel_a_mes)
dev.off()

# -----------------------------------------------------------------------------
# 3. PANEL B - Confronto composizione tumore vs linee 2D (deconvoluzione)
# -----------------------------------------------------------------------------
message("\n--- 3. Panel B: confronto composizione tumore vs linee ---")

meta_temp <- read_csv(SAMPLE_METADATA, show_col_types = FALSE)
sample_order <- meta_temp %>%
  arrange(patient, condition) %>%
  pull(sample_id)
bulk_prop$sample_id <- factor(bulk_prop$sample_id, levels = sample_order)
bulk_prop$condition <- factor(bulk_prop$condition, levels = c("tumor", "2D"))

panel_b_stack <- ggplot(
  bulk_prop,
  aes(x = sample_id, y = proportion, fill = cell_type)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = cell_colors) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  labs(
    title = "B \u2014 Deconvolution \u2014 Cell type composition",
    x = NULL, y = "Proportion", fill = "Cell type"
  ) +
  THEME_PGL +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text  = element_text(face = "bold", size = 11)
  )

# Boxplot mesenchimale paired
panel_b_mes <- ggplot(
  bulk_mes_prop,
  aes(x = condition, y = prop_mesenchymal * 100,
      color = condition, fill = condition)
) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(aes(group = patient), color = "grey40",   # paired lines
            alpha = 0.5, linewidth = 0.7) +
  scale_color_manual(values = COLORS_CONDITION) +
  scale_fill_manual(values  = COLORS_CONDITION) +
  labs(
    title = "B2 \u2014 Mesenchymal proportion",
    x = NULL, y = "% mesenchymal cells",
    color = NULL, fill = NULL
  ) +
  THEME_PGL + theme(legend.position = "none")

pdf(file.path(RESULTS_FIGURES, "panelB_deconvolution_comparison.pdf"),
    width = 18, height = 7)
print(panel_b_stack + panel_b_mes + plot_layout(widths = c(3, 1)))
dev.off()

# -----------------------------------------------------------------------------
# 4. PANEL C - Volcano con marker mesenchimali evidenziati
# -----------------------------------------------------------------------------
message("\n--- 4. Panel C: volcano DEG ---")

degs_all_plot <- degs_all %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < DEG_PADJ_THRESHOLD & log2FoldChange >  DEG_LFC_THRESHOLD ~ "Up in 2D",
      padj < DEG_PADJ_THRESHOLD & log2FoldChange < -DEG_LFC_THRESHOLD ~ "Up in tumor",
      TRUE ~ "NS"
    ),
    highlight = case_when(
      gene_symbol %in% MESENCHYMAL_MARKERS    ~ "Mesenchymal",
      gene_symbol %in% NEUROENDOCRINE_MARKERS ~ "Neuroendocrine",
      TRUE ~ NA_character_
    ),
    label = ifelse(!is.na(highlight), gene_symbol, NA_character_)
  )

panel_c <- ggplot(degs_all_plot,
                  aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = filter(degs_all_plot, is.na(highlight)),
             aes(color = sig), size = 0.6, alpha = 0.5) +
  geom_point(data = filter(degs_all_plot, !is.na(highlight)),
             aes(color = highlight), size = 2.5, alpha = 1) +
  geom_text_repel(
    data = filter(degs_all_plot, !is.na(label)),
    aes(label = label, color = highlight),
    size = 3.2, fontface = "italic",
    max.overlaps = 40, segment.color = "grey40"
  ) +
  scale_color_manual(
    values = c(
      "Up in 2D"      = "#4DBBD5",
      "Up in tumor"   = "#E64B35",
      "NS"            = "grey75",
      "Mesenchymal"   = "#FF7F0E",
      "Neuroendocrine" = "#9467BD"
    ),
    breaks = c("Up in 2D", "Up in tumor", "Mesenchymal", "Neuroendocrine")
  ) +
  geom_vline(xintercept = c(-DEG_LFC_THRESHOLD, DEG_LFC_THRESHOLD),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(DEG_PADJ_THRESHOLD),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  labs(
    title    = "C \u2014 Volcano plot",
    subtitle = sprintf("%d DEGs | Mesenchymal & neuroendocrine markers highlighted",
                       nrow(degs_filtered)),
    x = "log2 Fold Change (2D / tumor)",
    y = "-log10 (padj)",
    color = NULL
  ) +
  THEME_PGL

ggsave(file.path(RESULTS_FIGURES, "panelC_volcano_annotated.pdf"),
       panel_c, width = 10, height = 8)

# -----------------------------------------------------------------------------
# 5. PANEL D - Overlap DEG upregolati nelle 2D vs marker mesenchimali scRNA
# -----------------------------------------------------------------------------
message("\n--- 5. Panel D: overlap DEG vs marker mesenchimali scRNA ---")

# Geni upregolati nelle linee 2D
deg_up_2D <- degs_filtered %>%
  filter(log2FoldChange > DEG_LFC_THRESHOLD) %>%
  pull(gene_symbol) %>%
  na.omit() %>%
  unique()

# Geni downregolati nelle linee 2D (= upregolati nel tumore)
deg_up_tumor <- degs_filtered %>%
  filter(log2FoldChange < -DEG_LFC_THRESHOLD) %>%
  pull(gene_symbol) %>%
  na.omit() %>%
  unique()

# Marker mesenchimali dall'scRNA
# Identifica il nome della colonna cluster (puo' variare)
cluster_col <- intersect(c("cluster", "cell_type", "Cluster"),
                          colnames(scrna_markers))[1]
gene_col     <- intersect(c("gene", "gene_symbol", "Gene"),
                          colnames(scrna_markers))[1]

mesenchymal_pattern <- "(?i)(mesench|fibroblast|stromal|sustentacular|MSC)"
scrna_mes_markers <- scrna_markers %>%
  filter(grepl(mesenchymal_pattern, .data[[cluster_col]], perl = TRUE)) %>%
  pull(.data[[gene_col]]) %>%
  unique()

cat(sprintf("DEG up in 2D: %d\n", length(deg_up_2D)))
cat(sprintf("DEG up in tumor: %d\n", length(deg_up_tumor)))
cat(sprintf("Marker mesenchimali scRNA: %d\n", length(scrna_mes_markers)))

# Overlap
overlap_mes_2D    <- intersect(deg_up_2D,    scrna_mes_markers)
overlap_mes_tumor <- intersect(deg_up_tumor, scrna_mes_markers)

cat(sprintf("\nMarker mesenchimali upregolati nelle LINEE 2D: %d\n",
            length(overlap_mes_2D)))
cat(sprintf("Marker mesenchimali upregolati nel TUMORE: %d\n",
            length(overlap_mes_tumor)))
print(overlap_mes_2D)
print(overlap_mes_tumor)

write_csv(
  data.frame(gene = overlap_mes_2D),
  file.path(RESULTS_FIGURES, "overlap_mes_markers_upregulated_in_2D.csv")
)
write_csv(
  data.frame(gene = overlap_mes_tumor),
  file.path(RESULTS_FIGURES, "overlap_mes_markers_upregulated_in_tumor.csv")
)

# Venn diagram
venn_list <- list(
  "DEG up in 2D"         = deg_up_2D,
  "Marker mes. (scRNA)"  = scrna_mes_markers
)

pdf(file.path(RESULTS_FIGURES, "panelD_venn_deg_vs_mes_markers.pdf"),
    width = 7, height = 6)
ggvenn::ggvenn(
  venn_list,
  fill_color    = c("#4DBBD5", "#F39B7F"),
  stroke_size   = 0.5,
  text_size     = 5,
  set_name_size = 4
) +
  labs(title = "D \u2014 DEG vs mesenchymal marker overlap") +
  THEME_PGL
dev.off()

# Tabella geni in overlap (DEG up 2D + mesenchimali scRNA)
if (length(overlap_mes_2D) > 0) {
  overlap_table <- degs_filtered %>%
    filter(gene_symbol %in% overlap_mes_2D) %>%
    select(gene_symbol, log2FoldChange, padj, baseMean) %>%
    arrange(desc(log2FoldChange))

  pdf(file.path(RESULTS_FIGURES, "panelD_overlap_mes_genes_barplot.pdf"),
      width = 8, height = max(4, nrow(overlap_table) * 0.35 + 2))

  p_overlap <- ggplot(
    overlap_table,
    aes(x = reorder(gene_symbol, log2FoldChange), y = log2FoldChange,
        fill = log2FoldChange > 0)
  ) +
    geom_col() +
    scale_fill_manual(values = c("FALSE" = "#E64B35", "TRUE" = "#4DBBD5"),
                      guide = "none") +
    coord_flip() +
    labs(
      title = "Mesenchymal markers among DEGs",
      x = NULL, y = "log2FC (2D / tumor)"
    ) +
    THEME_PGL +
    theme(axis.text.y = element_text(face = "italic"))

  print(p_overlap)
  dev.off()
}

# -----------------------------------------------------------------------------
# 6. FIGURA RIASSUNTIVA - "le linee non sono 100% tumore"
# -----------------------------------------------------------------------------
message("\n--- 6. Figura riassuntiva ---")

# Proporzione "non-tumor" per ogni campione dalla deconvoluzione
# (tutto cio' che non e' classificato come neuroendocrino/chromaffin)
ne_pattern <- "(?i)(neuroend|chromaffin|chief|tumor|malignant|PGL|paraganglioma)"

nontumorprop_df <- bulk_prop %>%
  mutate(is_tumor_cell = grepl(ne_pattern, cell_type, perl = TRUE)) %>%
  group_by(sample_id, condition, patient) %>%
  summarise(
    prop_tumor    = sum(proportion[is_tumor_cell]),
    prop_nontumor = sum(proportion[!is_tumor_cell]),
    .groups = "drop"
  )

p_nontumor <- ggplot(
  nontumorprop_df,
  aes(x = reorder(sample_id, prop_nontumor),
      y = prop_nontumor * 100,
      fill = condition)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = COLORS_CONDITION) +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "grey20", linewidth = 0.7) +
  coord_flip() +
  labs(
    title    = "Non-tumor cell proportion per sample",
    subtitle = NULL,
    x = NULL, y = "% non-tumor cells (MuSiC estimate)",
    fill = "Condition"
  ) +
  THEME_PGL

ggsave(file.path(RESULTS_FIGURES, "summary_nontumor_proportion.pdf"),
       p_nontumor, width = 10, height = 7)

# -----------------------------------------------------------------------------
# 7. TABELLA RIASSUNTIVA FINALE
# -----------------------------------------------------------------------------
message("\n--- 7. Tabella riassuntiva ---")

summary_final <- data.frame(
  Analysis = c(
    "Total DEGs (2D vs tumor)",
    "DEGs up in 2D cell lines",
    "DEGs down in 2D lines (up in tumor)",
    "Mesenchymal scRNA markers among DEGs up in 2D",
    "Mesenchymal scRNA markers among DEGs up in tumor",
    "Samples with mesenchymal cells > 5% (tumor scRNA-seq)",
    "2D samples with non-tumor cells > 10% (deconvolution)"
  ),
  Valore = c(
    nrow(degs_filtered),
    length(deg_up_2D),
    length(deg_up_tumor),
    length(overlap_mes_2D),
    length(overlap_mes_tumor),
    sum(scrna_mes_prop$proportion_mesenchymal > 0.05),
    sum(nontumorprop_df$prop_nontumor > 0.1 &
          nontumorprop_df$condition == "2D")
  )
)

print(summary_final)
write_csv(summary_final, file.path(RESULTS_FIGURES, "summary_table.csv"))

message("\n=== 05_integration_figures.R completato ===")
message("Output in: ", RESULTS_FIGURES)
message("\n--- PIPELINE COMPLETATA ---")
