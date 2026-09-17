# =============================================================================
# 01_scrna_reference.R
# Caratterizzazione tipi cellulari nel tumore (scRNA-seq)
# + costruzione della matrice di riferimento per la deconvoluzione bulk (MuSiC)
#
# Annotazione: usa le annotazioni finali di Pasquale (Idents dell'oggetto Seurat)
#
# Input:  pasquale_s_obj_multi_finalAnnotation.rds
# Output: UMAP cell types, proporzioni mesenchimali, reference matrix per MuSiC
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(MuSiC)
  library(Biobase)
})

check_inputs(SCRNA_INTEGRATED_RDS)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO OGGETTO SEURAT DI PASQUALE
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento oggetto Seurat (Pasquale) ---")

seu <- readRDS(SCRNA_INTEGRATED_RDS)
cat(sprintf("Cellule: %d | Geni: %d\n", ncol(seu), nrow(seu)))

cat("\nColonne metadata disponibili:\n")
print(colnames(seu@meta.data))

cat("\nIdents (prime 20):\n")
print(head(Idents(seu), 20))

cat("\nValori unici Idents:\n")
print(sort(table(Idents(seu)), decreasing = TRUE))

# -----------------------------------------------------------------------------
# 2. ANNOTAZIONE: usa Idents di Pasquale come tipo cellulare
# -----------------------------------------------------------------------------
message("\n--- 2. Annotazione da Pasquale ---")

seu$celltype_pasquale <- as.character(Idents(seu))
CELLTYPE_COL  <- "celltype_pasquale"
SAMPLEID_COL  <- "orig.ident"

if (!SAMPLEID_COL %in% colnames(seu@meta.data)) {
  # Prova colonne alternative comuni
  alt_cols <- c("sample", "Sample", "orig.ident", "sampleID", "patient")
  found <- alt_cols[alt_cols %in% colnames(seu@meta.data)]
  if (length(found) == 0) stop("Nessuna colonna campione trovata nel metadata.")
  SAMPLEID_COL <- found[1]
  message(sprintf("Usata colonna campione: %s", SAMPLEID_COL))
}

cat(sprintf("\nColonna campione: %s\n", SAMPLEID_COL))
cat("Campioni:\n")
print(table(seu@meta.data[[SAMPLEID_COL]]))

cat("\nTipi cellulari (Idents di Pasquale):\n")
print(sort(table(seu$celltype_pasquale), decreasing = TRUE))

cat("\nAnnotazione per campione:\n")
print(table(seu@meta.data[[SAMPLEID_COL]], seu$celltype_pasquale))

# Colori per i tipi cellulari trovati
cell_types_found <- sort(unique(seu$celltype_pasquale))
n_types <- length(cell_types_found)
COLORS_CELLTYPE <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(min(12, n_types), "Paired"))(n_types),
  cell_types_found
)

# -----------------------------------------------------------------------------
# 3. UMAP - VISUALIZZAZIONE TIPI CELLULARI
# -----------------------------------------------------------------------------
message("\n--- 3. UMAP tipi cellulari ---")

p_umap_celltypes <- DimPlot(
  seu, reduction = "umap", group.by = CELLTYPE_COL,
  label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.3,
  cols = COLORS_CELLTYPE
) +
  labs(title = "Cell types (Pasquale annotation)") +
  THEME_PGL

p_umap_sample <- DimPlot(
  seu, reduction = "umap", group.by = SAMPLEID_COL,
  label = FALSE, pt.size = 0.3
) +
  labs(title = "Sample of origin") +
  THEME_PGL

pdf(file.path(RESULTS_SCRNA, "umap_celltypes_overview.pdf"), width = 14, height = 6)
print(p_umap_celltypes | p_umap_sample)
dev.off()

# -----------------------------------------------------------------------------
# 4. PROPORZIONI TIPI CELLULARI PER CAMPIONE
# -----------------------------------------------------------------------------
message("\n--- 4. Proporzioni cellulari per campione ---")

prop_df <- seu@meta.data %>%
  as.data.frame() %>%
  group_by(.data[[SAMPLEID_COL]], .data[[CELLTYPE_COL]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[SAMPLEID_COL]]) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
colnames(prop_df)[1:2] <- c("sample", "cell_type")

write_csv(prop_df, file.path(RESULTS_SCRNA, "celltype_proportions_per_sample.csv"))

# Identifica cluster mesenchimali/stromali
mesenchymal_pattern <- "(?i)(mesench|fibroblast|stromal|sustentacular|MSC)"
prop_df$is_mesenchymal <- grepl(mesenchymal_pattern, prop_df$cell_type, perl = TRUE)

prop_mes <- prop_df %>%
  group_by(sample) %>%
  summarise(
    proportion_mesenchymal = sum(proportion[is_mesenchymal]),
    n_mesenchymal          = sum(n[is_mesenchymal]),
    .groups = "drop"
  ) %>%
  arrange(desc(proportion_mesenchymal))

cat("\nProporzione cellule mesenchimali/stromali per campione:\n")
print(prop_mes)
write_csv(prop_mes, file.path(RESULTS_SCRNA, "mesenchymal_proportions.csv"))

# Barplot composizione cellulare per campione
p_prop <- ggplot(prop_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = COLORS_CELLTYPE) +
  labs(
    title = "Cell type composition per sample (tumor scRNA-seq)",
    x = "Sample", y = "Proportion", fill = "Cell type"
  ) +
  THEME_PGL +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(RESULTS_SCRNA, "celltype_proportions_barplot.pdf"),
       p_prop, width = 12, height = 7)

# Barplot % mesenchimali per campione
p_mes <- ggplot(prop_mes,
                aes(x = reorder(sample, proportion_mesenchymal),
                    y = proportion_mesenchymal * 100)) +
  geom_col(fill = "#4DBBD5") +
  geom_text(aes(label = sprintf("%.1f%%", proportion_mesenchymal * 100)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "% mesenchymal/stromal cells per sample (tumor scRNA-seq)",
    x = "Sample", y = "% mesenchymal cells"
  ) +
  ylim(0, max(prop_mes$proportion_mesenchymal * 100, na.rm = TRUE) * 1.2) +
  THEME_PGL

ggsave(file.path(RESULTS_SCRNA, "mesenchymal_proportion_per_sample.pdf"),
       p_mes, width = 8, height = 6)

# -----------------------------------------------------------------------------
# 5. MARKER DEI TIPI CELLULARI
# -----------------------------------------------------------------------------
message("\n--- 5. Marker per tipo cellulare ---")

if (file.exists(MARKER_FILE)) {
  message("Trovato file marker di Pasquale - uso quelli esistenti.")
  markers <- readxl::read_xlsx(MARKER_FILE)
  write_csv(markers, file.path(RESULTS_SCRNA, "markers_per_celltype.csv"))
} else {
  message("Calcolo marker de novo con FindAllMarkers...")
  Idents(seu) <- CELLTYPE_COL
  markers <- FindAllMarkers(
    seu,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.5,
    test.use        = "wilcox"
  )
  write_csv(markers, file.path(RESULTS_SCRNA, "markers_per_celltype.csv"))
}

# -----------------------------------------------------------------------------
# 6. COSTRUZIONE REFERENCE MATRIX PER MuSiC
# -----------------------------------------------------------------------------
message("\n--- 6. Costruzione ExpressionSet per MuSiC ---")

tryCatch({
  counts_mat <- GetAssayData(seu, assay = "RNA", slot = "counts")
}, error = function(e) {
  counts_mat <<- GetAssayData(seu, assay = "RNA", layer = "counts")
})

cat(sprintf("Dimensioni count matrix: %d geni x %d cellule\n",
            nrow(counts_mat), ncol(counts_mat)))

pdata <- data.frame(
  row.names    = colnames(counts_mat),
  cellType     = seu@meta.data[[CELLTYPE_COL]],
  sampleID     = seu@meta.data[[SAMPLEID_COL]],
  stringsAsFactors = FALSE
)

# Rimuovi cellule senza annotazione valida
valid_cells <- !is.na(pdata$cellType) & pdata$cellType != "" &
               !pdata$cellType %in% c("unknown", "Unknown", "NA")
pdata      <- pdata[valid_cells, ]
counts_mat <- counts_mat[, valid_cells]

cat(sprintf("Cellule con annotazione valida: %d\n", ncol(counts_mat)))
cat("Tipi cellulari nel reference MuSiC:\n")
print(sort(table(pdata$cellType), decreasing = TRUE))

scrna_eset <- ExpressionSet(
  assayData = as.matrix(counts_mat),
  phenoData = new("AnnotatedDataFrame", data = pdata)
)

saveRDS(scrna_eset,
        file.path(RESULTS_SCRNA, "scrna_expressionset_for_music.RDS"))
message("ExpressionSet salvato per MuSiC.")

# -----------------------------------------------------------------------------
# 7. FEATUREPLOT MARKER MESENCHIMALI
# -----------------------------------------------------------------------------
message("\n--- 7. FeaturePlot marker mesenchimali ---")

mes_available <- MESENCHYMAL_MARKERS[MESENCHYMAL_MARKERS %in% rownames(seu)]
cat(sprintf("Marker mesenchimali disponibili: %d/%d\n",
            length(mes_available), length(MESENCHYMAL_MARKERS)))

if (length(mes_available) > 0) {
  plot_genes <- head(mes_available, 12)
  pdf(file.path(RESULTS_SCRNA, "featureplot_mesenchymal_markers.pdf"),
      width = 16, height = 12)
  print(
    FeaturePlot(seu, features = plot_genes,
                reduction = "umap", pt.size = 0.2,
                ncol = 4, order = TRUE) &
      theme(legend.position = "right", text = element_text(size = 9))
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# 8. SALVATAGGIO METADATI ESTRATTI
# -----------------------------------------------------------------------------
message("\n--- 8. Salvataggio metadati ---")

metadata_export <- seu@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_barcode") %>%
  select(cell_barcode,
         sample    = all_of(SAMPLEID_COL),
         cell_type = all_of(CELLTYPE_COL),
         any_of(c("nCount_RNA", "nFeature_RNA",
                  "percent.mt", "seurat_clusters")))

write_csv(metadata_export,
          file.path(RESULTS_SCRNA, "scrna_metadata_export.csv"))

message("\n=== 01_scrna_reference.R completato ===")
message("Output in: ", RESULTS_SCRNA)
