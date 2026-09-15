# =============================================================================
# 00_config.R
# Configurazione globale della pipeline tumore vs linee primarie 2D
# Sourcato da tutti gli altri script con: source("00_config.R")
# =============================================================================

# -----------------------------------------------------------------------------
# PERCORSI BASE
# Rileva automaticamente se gira su cluster (Terastat) o in locale
# -----------------------------------------------------------------------------
.node <- Sys.info()[["nodename"]]
.on_cluster <- grepl("terastat|login|compute", .node, ignore.case = TRUE)

if (.on_cluster) {
  # --- CLUSTER (Terastat, Sapienza) ---
  BASE_DIR     <- "/lustre/home/gfiscon/projects/PGL"
  PIPELINE_DIR <- file.path(BASE_DIR, "pipeline")
  DATA_DIR     <- file.path(PIPELINE_DIR, "data")
  RESULTS_DIR  <- file.path(BASE_DIR, "results")
  N_CORES      <- 8   # adatta in base alle CPU richieste nel job SLURM

  SCRNA_INTEGRATED_RDS <- file.path(BASE_DIR,
    "data/scrna_integrated/PGL_PFE_3_integrated_x_DE_complete_annotation.Rds")
  SCEVAN_DIR      <- file.path(BASE_DIR, "data/scevan")
  BULK_SALMON_DIR <- file.path(BASE_DIR, "data/bulk_salmon")
  MARKER_FILE     <- file.path(BASE_DIR, "reference/marker_per_cluster.xlsx")

} else {
  # --- LOCALE (Mac) ---
  BASE_DIR <- "/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/PGL"
  PIPELINE_DIR <- file.path(BASE_DIR, "analisi_fede/pipeline_linee_vs_tumore")
  DATA_DIR     <- file.path(PIPELINE_DIR, "data")
  RESULTS_DIR  <- file.path(PIPELINE_DIR, "results")
  N_CORES      <- 4

  SCRNA_INTEGRATED_RDS <- file.path(BASE_DIR,
    "analisi_alessio/PGL_PFE_3_integrated_x_DE_complete_annotation.Rds")
  SCEVAN_DIR      <- file.path(BASE_DIR, "analisi_alessio/scevan_iterato_tutti_campioni")
  BULK_SALMON_DIR <- file.path(BASE_DIR, "analisi_pasquale/bulk RNA/RNA counts")
  MARKER_FILE     <- file.path(BASE_DIR,
    "analisi_pasquale/scRNA PGL/results/marker_per_cluster.xlsx")
}

message(sprintf("Ambiente: %s | BASE_DIR: %s", ifelse(.on_cluster, "CLUSTER", "LOCALE"), BASE_DIR))

SCEVAN_SAMPLES <- c("PC190", "PTJ173", "PTJ184", "PTJ185",
                    "PV158 BIS", "PV180", "PV181", "PV193")

# Input - Metadata campioni
SAMPLE_METADATA <- file.path(DATA_DIR, "sample_metadata.csv")

# -----------------------------------------------------------------------------
# CARTELLE OUTPUT (create automaticamente se non esistono)
# -----------------------------------------------------------------------------
RESULTS_SCRNA    <- file.path(RESULTS_DIR, "scrna_reference")
RESULTS_BULK     <- file.path(RESULTS_DIR, "bulk_preprocessing")
RESULTS_DEGS     <- file.path(RESULTS_DIR, "bulk_degs")
RESULTS_DEGS_FIG <- file.path(RESULTS_DIR, "bulk_degs/figures")
RESULTS_DECONV   <- file.path(RESULTS_DIR, "deconvolution")
RESULTS_DECONV_FIG <- file.path(RESULTS_DIR, "deconvolution/figures")
RESULTS_FIGURES  <- file.path(RESULTS_DIR, "figures")

for (d in c(RESULTS_SCRNA, RESULTS_BULK, RESULTS_DEGS, RESULTS_DEGS_FIG,
            RESULTS_DECONV, RESULTS_DECONV_FIG, RESULTS_FIGURES)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# PARAMETRI ANALISI
# -----------------------------------------------------------------------------
# DESeq2
DEG_PADJ_THRESHOLD <- 0.05
DEG_LFC_THRESHOLD  <- 1.0        # log2FC assoluto
DEG_MIN_COUNTS     <- 10         # read minime per gene
DEG_MIN_SAMPLES    <- 3          # in almeno N campioni

# Deconvoluzione
DECONV_N_MARKERS   <- 200        # marker per tipo cellulare per MuSiC

# Seed riproducibilità
SEED <- 42
set.seed(SEED)

# -----------------------------------------------------------------------------
# MARKER MESENCHIMALI DI RIFERIMENTO (letteratura)
# Usati per check espliciti nei DEG e nelle firme scRNA
# -----------------------------------------------------------------------------
MESENCHYMAL_MARKERS <- c(
  "VIM", "FN1", "S100A4", "ACTA2", "FAP",
  "THY1", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2",
  "COL3A1", "POSTN", "TWIST1", "SNAI1", "SNAI2",
  "ZEB1", "ZEB2", "CDH2", "MMP2", "MMP9"
)

# Marker neuroendocrini/tumorali PGL (cellule cromaffini)
NEUROENDOCRINE_MARKERS <- c(
  "CHGA", "CHGB", "SYP", "TH", "DBH",
  "PHOX2B", "PHOX2A", "RET", "HAND2",
  "PNMT", "NPY", "SCG2", "SCG3"
)

# Marker sustentacolari (cellule di supporto specifiche PGL)
SUSTENTACULAR_MARKERS <- c(
  "S100B", "SOX10", "GFAP", "NGFR", "NES",
  "PLP1", "CDH19", "ERBB3", "L1CAM"
)

# Marker endoteliali
ENDOTHELIAL_MARKERS <- c(
  "PECAM1", "VWF", "CDH5", "CLDN5", "ENG",
  "ESAM", "PTPRB", "ROBO4", "KDR"
)

# Marker immunitari (macrofagi, linfociti T, linfociti B)
IMMUNE_MARKERS <- c(
  "PTPRC", "CD68", "C1QA", "C1QB", "MRC1",
  "CD3D", "CD3E", "CD8A", "CD4",
  "CD79A", "MS4A1", "NKG7"
)

# -----------------------------------------------------------------------------
# TEMA GGPLOT CONDIVISO
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))

THEME_PGL <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey95"),
    legend.position   = "right",
    plot.title        = element_text(face = "bold", size = 13)
  )

# Palette condizionale
COLORS_CONDITION <- c("tumor" = "#E64B35", "2D" = "#4DBBD5")

# -----------------------------------------------------------------------------
# FUNZIONE DI CHECK INPUT
# Verifica che i file richiesti esistano prima di partire
# -----------------------------------------------------------------------------
check_inputs <- function(...) {
  files <- c(...)
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("File di input mancanti:\n", paste(" -", missing, collapse = "\n"))
  }
  message("OK: tutti i file di input trovati.")
}

message("00_config.R caricato.")
