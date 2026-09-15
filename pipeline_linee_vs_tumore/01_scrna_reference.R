# =============================================================================
# 01_scrna_reference.R
# Caratterizzazione tipi cellulari nel tumore (scRNA-seq)
# + costruzione della matrice di riferimento per la deconvoluzione bulk (MuSiC)
#
# Annotazione: marker-based con AddModuleScore() usando i marker definiti in
# 00_config.R (NEUROENDOCRINE_MARKERS, MESENCHYMAL_MARKERS).
# Non richiede SingleR o celldex.
#
# Input:  PGL_PFE_3_integrated_x_DE_complete_annotation.Rds (oggetto Seurat)
#         scevan_iterato_tutti_campioni/*/scevan_PGL.RDS (label tumor/normal)
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

SAMPLEID_COL <- "orig.ident"   # colonna campione (quasi certamente corretta)

check_inputs(SCRNA_INTEGRATED_RDS)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO OGGETTO SEURAT INTEGRATO
# Nota: file da 7.9 GB - richiede ~16-20 GB di RAM
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento oggetto Seurat integrato ---")
message("Questo passo richiede diversi minuti e ~16-20 GB di RAM...")

seu <- readRDS(SCRNA_INTEGRATED_RDS)
cat(sprintf("Cellule: %d | Geni: %d\n", ncol(seu), nrow(seu)))

cat("\nColonne metadata disponibili:\n")
print(colnames(seu@meta.data))

if (!SAMPLEID_COL %in% colnames(seu@meta.data)) {
  stop(sprintf("Colonna '%s' non trovata. Colonne disponibili:\n%s",
               SAMPLEID_COL, paste(colnames(seu@meta.data), collapse = ", ")))
}
cat("\nCampioni:\n")
print(table(seu@meta.data[[SAMPLEID_COL]]))

# -----------------------------------------------------------------------------
# 2. ANNOTAZIONE MARKER-BASED (Neuroendocrine / Mesenchymal / Other)
# Usa AddModuleScore con i marker definiti in 00_config.R.
# Piu' appropriato per PGL rispetto a SingleR con riferimenti generici.
# -----------------------------------------------------------------------------
message("\n--- 2. Annotazione marker-based ---")

# Filtra marker presenti nel dataset
marker_lists <- list(
  Neuroendocrine  = NEUROENDOCRINE_MARKERS,
  Mesenchymal     = MESENCHYMAL_MARKERS,
  Sustentacular   = SUSTENTACULAR_MARKERS,
  Endothelial     = ENDOTHELIAL_MARKERS,
  Immune          = IMMUNE_MARKERS
)

avail_lists <- lapply(marker_lists, function(m) m[m %in% rownames(seu)])

for (nm in names(avail_lists)) {
  cat(sprintf("%-16s: %d/%d marker disponibili\n",
              nm, length(avail_lists[[nm]]), length(marker_lists[[nm]])))
}

# Rimuovi popolazioni con meno di 3 marker disponibili
avail_lists <- Filter(function(m) length(m) >= 3, avail_lists)
if (length(avail_lists) == 0) stop("Nessuna popolazione con abbastanza marker.")

# Calcola module scores per ogni popolazione
score_names <- paste0(names(avail_lists), "_score1")
for (nm in names(avail_lists)) {
  seu <- AddModuleScore(seu,
                        features = list(avail_lists[[nm]]),
                        name     = paste0(nm, "_score"),
                        assay    = "RNA")
}
# AddModuleScore aggiunge suffisso "1": Neuroendocrine_score1, etc.

# Assegna il tipo cellulare: vince il punteggio piu' alto (purche' > 0)
score_mat <- as.data.frame(seu@meta.data[, score_names, drop = FALSE])
max_score <- apply(score_mat, 1, max)
max_pop   <- names(avail_lists)[apply(score_mat, 1, which.max)]

seu$marker_celltype <- ifelse(max_score > 0, max_pop, "Other")

CELLTYPE_COL <- "marker_celltype"

cat("\nDistribuzione annotazione marker-based:\n")
print(sort(table(seu$marker_celltype), decreasing = TRUE))

cat("\nAnnotazione per campione:\n")
print(table(seu@meta.data[[SAMPLEID_COL]], seu$marker_celltype))

# Palette colori per le 6 popolazioni
COLORS_CELLTYPE <- c(
  "Neuroendocrine" = "#E64B35",
  "Mesenchymal"    = "#4DBBD5",
  "Sustentacular"  = "#00A087",
  "Endothelial"    = "#F39B7F",
  "Immune"         = "#8491B4",
  "Other"          = "#B09C85"
)

# -----------------------------------------------------------------------------
# 3. AGGIUNTA LABEL SCEVAN (tumor / normal) DAI FILE PER-SAMPLE
# I per-sample scevan_PGL.RDS contengono: class (tumor/normal), subclone
# -----------------------------------------------------------------------------
message("\n--- 2. Aggiunta label SCEVAN (tumor/normal) ---")

scevan_labels <- lapply(SCEVAN_SAMPLES, function(sample_name) {
  scev_path <- file.path(SCEVAN_DIR, sample_name, "scevan_PGL.RDS")
  if (!file.exists(scev_path)) {
    warning(sprintf("SCEVAN file non trovato: %s", scev_path))
    return(NULL)
  }
  df <- readRDS(scev_path)  # data.frame: class, confidentNormal, subclone
  df$cell_barcode <- rownames(df)
  df$sample       <- sample_name
  return(df)
})

scevan_df <- bind_rows(Filter(Negate(is.null), scevan_labels))
cat(sprintf("Barcode con label SCEVAN: %d\n", nrow(scevan_df)))
cat("Distribuzione classe SCEVAN:\n")
print(table(scevan_df$class))

# Aggiungi la label SCEVAN al metadata del Seurat
# Nota: i barcode nel Seurat integrato possono avere suffissi (es. _1, _2)
# Il matching viene fatto con un approccio flessibile
seu$scevan_class <- NA_character_
matched <- match(rownames(seu@meta.data), scevan_df$cell_barcode)
seu$scevan_class[!is.na(matched)] <- scevan_df$class[matched[!is.na(matched)]]

cat(sprintf("Cellule con label SCEVAN assegnata: %d / %d\n",
            sum(!is.na(seu$scevan_class)), ncol(seu)))

# -----------------------------------------------------------------------------
# 4. UMAP - VISUALIZZAZIONE TIPI CELLULARI
# -----------------------------------------------------------------------------
message("\n--- 4. UMAP tipi cellulari ---")

Idents(seu) <- CELLTYPE_COL

p_umap_celltypes <- DimPlot(
  seu, reduction = "umap", group.by = CELLTYPE_COL,
  label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.3,
  cols = COLORS_CELLTYPE
) +
  labs(title = "Cell types (marker-based)") +
  THEME_PGL

p_umap_sample <- DimPlot(
  seu, reduction = "umap", group.by = SAMPLEID_COL,
  label = FALSE, pt.size = 0.3
) +
  labs(title = "Sample of origin") +
  THEME_PGL

p_umap_scevan <- DimPlot(
  seu, reduction = "umap", group.by = "scevan_class",
  label = FALSE, pt.size = 0.3, na.value = "grey90",
  cols = c("tumor" = "#E64B35", "normal" = "#4DBBD5")
) +
  labs(title = "SCEVAN: tumor vs normal") +
  THEME_PGL

pdf(file.path(RESULTS_SCRNA, "umap_celltypes_overview.pdf"), width = 18, height = 6)
print(p_umap_celltypes | p_umap_sample | p_umap_scevan)
dev.off()

# UMAP con tutti i module scores come feature continua
score_plots <- lapply(names(avail_lists), function(nm) {
  FeaturePlot(seu, features = paste0(nm, "_score1"),
              reduction = "umap", pt.size = 0.2, order = TRUE) +
    labs(title = paste(nm, "score")) + THEME_PGL
})

pdf(file.path(RESULTS_SCRNA, "umap_module_scores.pdf"),
    width = 7 * min(3, length(score_plots)),
    height = 7 * ceiling(length(score_plots) / 3))
print(patchwork::wrap_plots(score_plots, ncol = 3))
dev.off()

# -----------------------------------------------------------------------------
# 5. PROPORZIONI TIPI CELLULARI PER CAMPIONE
# Focus su cellule mesenchimali nel tumore originale
# -----------------------------------------------------------------------------
message("\n--- 5. Proporzioni cellulari per campione ---")

prop_df <- seu@meta.data %>%
  as.data.frame() %>%
  group_by(.data[[SAMPLEID_COL]], .data[[CELLTYPE_COL]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(.data[[SAMPLEID_COL]]) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
colnames(prop_df)[1:2] <- c("sample", "cell_type")

write_csv(prop_df, file.path(RESULTS_SCRNA, "celltype_proportions_per_sample.csv"))

# Identifica cluster mesenchimali (cerca pattern nel nome del tipo cellulare)
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

cat("\nProporzione cellule mesenchimali per campione:\n")
print(prop_mes)
write_csv(prop_mes, file.path(RESULTS_SCRNA, "mesenchymal_proportions.csv"))

# Barplot: composizione cellulare per campione
n_types <- length(unique(prop_df$cell_type))
cell_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_types),
  unique(prop_df$cell_type)
)

p_prop <- ggplot(prop_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_colors) +
  labs(
    title = "Cell type composition per sample (tumor scRNA-seq)",
    x = "Sample", y = "Proportion", fill = "Cell type"
  ) +
  THEME_PGL +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(RESULTS_SCRNA, "celltype_proportions_barplot.pdf"),
       p_prop, width = 12, height = 7)

# Barplot specifico: % cellule mesenchimali per campione
p_mes <- ggplot(prop_mes,
                aes(x = reorder(sample, proportion_mesenchymal),
                    y = proportion_mesenchymal * 100)) +
  geom_col(fill = "#F39B7F") +
  geom_text(aes(label = sprintf("%.1f%%", proportion_mesenchymal * 100)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "% mesenchymal cells per sample (tumor scRNA-seq)",
    x = "Sample", y = "% mesenchymal cells"
  ) +
  ylim(0, max(prop_mes$proportion_mesenchymal * 100) * 1.2) +
  THEME_PGL

ggsave(file.path(RESULTS_SCRNA, "mesenchymal_proportion_per_sample.pdf"),
       p_mes, width = 8, height = 6)

# -----------------------------------------------------------------------------
# 6. MARKER DEI TIPI CELLULARI (o carica quelli gia' calcolati da Pasquale)
# -----------------------------------------------------------------------------
message("\n--- 6. Marker per tipo cellulare ---")

marker_file_pasquale <- MARKER_FILE

if (file.exists(marker_file_pasquale)) {
  # Usa i marker gia' calcolati da Pasquale
  message("Trovato file marker di Pasquale - uso quelli esistenti.")
  markers <- readxl::read_xlsx(marker_file_pasquale)
  write_csv(markers, file.path(RESULTS_SCRNA, "markers_per_celltype.csv"))
} else {
  # Calcola de novo (piu' lento)
  message("Calcolo marker de novo con FindAllMarkers...")
  Idents(seu) <- CELLTYPE_COL
  markers <- FindAllMarkers(
    seu,
    only.pos     = TRUE,
    min.pct      = 0.25,
    logfc.threshold = 0.5,
    test.use     = "wilcox"
  )
  write_csv(markers, file.path(RESULTS_SCRNA, "markers_per_celltype.csv"))
}

cat("Marker calcolati/caricati per i tipi cellulari:\n")
print(table(markers[[ifelse("cluster" %in% colnames(markers), "cluster", colnames(markers)[1])]]))

# -----------------------------------------------------------------------------
# 7. COSTRUZIONE REFERENCE MATRIX PER MuSiC
# Richiede: matrice counts (geni x cellule), label tipo cellulare, label paziente
# MuSiC usa la variabilita' inter-soggetto per pesare i geni nel reference
# -----------------------------------------------------------------------------
message("\n--- 7. Costruzione ExpressionSet per MuSiC ---")

# Estrai la count matrix raw (slot RNA, layer counts)
# In Seurat v5 la sintassi puo' essere diversa - gestito entrambi i casi
tryCatch({
  counts_mat <- GetAssayData(seu, assay = "RNA", slot = "counts")
}, error = function(e) {
  counts_mat <<- GetAssayData(seu, assay = "RNA", layer = "counts")
})

cat(sprintf("Dimensioni count matrix: %d geni x %d cellule\n",
            nrow(counts_mat), ncol(counts_mat)))

# Phenotype data per ExpressionSet MuSiC
pdata <- data.frame(
  row.names    = colnames(counts_mat),
  cellType     = seu@meta.data[[CELLTYPE_COL]],
  sampleID     = seu@meta.data[[SAMPLEID_COL]],
  stringsAsFactors = FALSE
)

# Rimuovi cellule senza annotazione
valid_cells <- !is.na(pdata$cellType) & pdata$cellType != "" &
               pdata$cellType != "unknown" & pdata$cellType != "Unknown"
pdata       <- pdata[valid_cells, ]
counts_mat  <- counts_mat[, valid_cells]

cat(sprintf("Cellule con annotazione valida: %d\n", ncol(counts_mat)))
cat("Tipi cellulari nel reference:\n")
print(sort(table(pdata$cellType), decreasing = TRUE))

# Costruisci ExpressionSet
pheno_data   <- new("AnnotatedDataFrame", data = pdata)
scrna_eset   <- ExpressionSet(
  assayData    = as.matrix(counts_mat),
  phenoData    = pheno_data
)

saveRDS(scrna_eset,
        file.path(RESULTS_SCRNA, "scrna_expressionset_for_music.RDS"))

message("ExpressionSet salvato per MuSiC.")

# -----------------------------------------------------------------------------
# 8. FEATURE PLOT MARKER MESENCHIMALI NEL SEURAT
# Visualizza dove si esprimono i marker mesenchimali nell'UMAP del tumore
# -----------------------------------------------------------------------------
message("\n--- 8. FeaturePlot marker mesenchimali ---")

# Controlla quali marker mesenchimali sono nel Seurat
mes_available <- MESENCHYMAL_MARKERS[
  MESENCHYMAL_MARKERS %in% rownames(seu)
]
cat(sprintf("Marker mesenchimali disponibili nel Seurat: %d/%d\n",
            length(mes_available), length(MESENCHYMAL_MARKERS)))

if (length(mes_available) > 0) {
  # Griglia di FeaturePlot per i primi 12 marker disponibili
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
# 9. SALVATAGGIO METADATI ESTRATTI (senza object pesante)
# -----------------------------------------------------------------------------
message("\n--- 9. Salvataggio metadati ---")

metadata_export <- seu@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_barcode") %>%
  select(cell_barcode,
         sample     = all_of(SAMPLEID_COL),
         cell_type  = all_of(CELLTYPE_COL),
         scevan_class,
         any_of(score_names),
         any_of(c("nCount_RNA", "nFeature_RNA",
                  "percent.mt", "seurat_clusters")))

write_csv(metadata_export,
          file.path(RESULTS_SCRNA, "scrna_metadata_export.csv"))

message("\n=== 01_scrna_reference.R completato ===")
message("Output in: ", RESULTS_SCRNA)
