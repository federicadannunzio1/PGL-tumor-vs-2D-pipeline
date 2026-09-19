# =============================================================================
# 06_broader_annotation.R
# Crea annotazione ampia sull'oggetto Seurat di Pasquale
# + costruzione nuovo ExpressionSet per MuSiC (senza TAM-like macrophages)
#
# Input:  pasquale_s_obj_multi_finalAnnotation.rds
# Output: ExpressionSet con broader_annotation per MuSiC
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Biobase)
  library(dplyr)
  library(readr)
})

check_inputs(SCRNA_INTEGRATED_RDS)

# Output directory
RESULTS_DECONV_BROAD     <- file.path(RESULTS_DIR, "deconvolution_broader")
RESULTS_DECONV_BROAD_FIG <- file.path(RESULTS_DIR, "deconvolution_broader/figures")
dir.create(RESULTS_DECONV_BROAD,     showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DECONV_BROAD_FIG, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO SEURAT
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento Seurat ---")

seu <- readRDS(SCRNA_INTEGRATED_RDS)
seu$celltype_pasquale <- as.character(Idents(seu))

cat("Tipi cellulari originali:\n")
print(sort(table(seu$celltype_pasquale), decreasing = TRUE))

# -----------------------------------------------------------------------------
# 2. MAPPING ANNOTAZIONE AMPIA (da tipi_cellulari_annotazione_pasquale.xlsx)
# -----------------------------------------------------------------------------
message("\n--- 2. Mapping annotazione ampia ---")

broader_map <- c(
  "Pericyte-like/mural cells"               = "Stroma",
  "Pericyte-like stromal-mesenchymal cells"  = "Stroma",
  "PLVAP+ permeable endothelial cells"       = "Endothelia",
  "PLVAP+ angiogenic endothelial cells"      = "Endothelia",
  "Arterial-like endothelial cells"          = "Endothelia",
  "Venous/venular endothelial cells"         = "Endothelia",
  "ACKR1+ venular endothelial cells"         = "Endothelia",
  "Inflammatory monocytes/macrophages"       = "Lymphocytes",
  "IL7R+ activated/memory T cells"           = "Lymphocytes",
  "Cycling myeloid cells"                    = "Lymphocytes",
  "TAM-like macrophages"                     = "Lymphocytes",
  "NK cells"                                 = "Lymphocytes",
  "Activated B cells"                        = "Lymphocytes",
  "Antigen-presenting myeloid cells"         = "Lymphocytes",
  "Smooth muscle-like cells"                 = "Smooth muscle-like cells",
  "Fibroblasts"                              = "Fibroblasts",
  "Myofibroblasts"                           = "Fibroblasts",
  "Neuroendocrine"                           = "Neuroendocrine",
  "Neuroendocrine PV158"                     = "Neuroendocrine",
  "Sustentacular"                            = "Sustentacular"
)

# Applica mapping
seu@meta.data$broader_annotation <- unname(broader_map[seu$celltype_pasquale])

# Controlla tipi non mappati
unmapped <- is.na(seu$broader_annotation)
if (any(unmapped)) {
  message("ATTENZIONE: tipi cellulari non mappati:")
  print(table(seu$celltype_pasquale[unmapped]))
  seu$broader_annotation[unmapped] <- "Other"
}

cat("\nAnnotazione ampia (tutte le cellule):\n")
print(sort(table(seu$broader_annotation), decreasing = TRUE))

# Salva mapping come CSV per riferimento
mapping_df <- data.frame(
  celltype_pasquale  = names(broader_map),
  broader_annotation = unname(broader_map),
  stringsAsFactors   = FALSE
)
write_csv(mapping_df, file.path(RESULTS_SCRNA, "broader_annotation_mapping.csv"))

# -----------------------------------------------------------------------------
# 3. RIMOZIONE TAM-like macrophages DAL REFERENCE
# Motivo: evitare falsi positivi nella deconvoluzione
# -----------------------------------------------------------------------------
message("\n--- 3. Rimozione TAM-like macrophages ---")

n_before <- ncol(seu)
seu_notam <- subset(seu, celltype_pasquale != "TAM-like macrophages")
n_after <- ncol(seu_notam)
cat(sprintf("Cellule rimosse (TAM): %d (da %d a %d)\n",
            n_before - n_after, n_before, n_after))

cat("\nAnnotazione ampia dopo rimozione TAM:\n")
print(sort(table(seu_notam$broader_annotation), decreasing = TRUE))

# -----------------------------------------------------------------------------
# 4. COSTRUZIONE ExpressionSet CON BROADER ANNOTATION
# -----------------------------------------------------------------------------
message("\n--- 4. Costruzione ExpressionSet broader per MuSiC ---")

SAMPLEID_COL <- "orig.ident"

tryCatch({
  counts_mat <- GetAssayData(seu_notam, assay = "RNA", slot = "counts")
}, error = function(e) {
  counts_mat <<- GetAssayData(seu_notam, assay = "RNA", layer = "counts")
})

cat(sprintf("Dimensioni count matrix: %d geni x %d cellule\n",
            nrow(counts_mat), ncol(counts_mat)))

pdata <- data.frame(
  row.names = colnames(counts_mat),
  cellType  = seu_notam$broader_annotation,
  sampleID  = seu_notam@meta.data[[SAMPLEID_COL]],
  stringsAsFactors = FALSE
)

# Rimuovi cellule senza annotazione valida
valid_cells <- !is.na(pdata$cellType) & pdata$cellType != "" &
               !pdata$cellType %in% c("unknown", "Unknown", "NA", "Other")
pdata      <- pdata[valid_cells, ]
counts_mat <- counts_mat[, valid_cells]

cat(sprintf("Cellule con annotazione valida: %d\n", ncol(counts_mat)))
cat("Tipi cellulari nel reference MuSiC (broader, no TAM):\n")
print(sort(table(pdata$cellType), decreasing = TRUE))

scrna_eset_broad <- ExpressionSet(
  assayData = as.matrix(counts_mat),
  phenoData = new("AnnotatedDataFrame", data = pdata)
)

saveRDS(scrna_eset_broad,
        file.path(RESULTS_SCRNA, "scrna_expressionset_broader_for_music.RDS"))

message("\n=== 06_broader_annotation.R completato ===")
message("ExpressionSet broader salvato in: ", RESULTS_SCRNA)
