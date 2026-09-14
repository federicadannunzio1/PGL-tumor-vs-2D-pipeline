# =============================================================================
# install_packages.R
# Installa tutti i pacchetti R necessari nella libreria personale su Terastat
#
# Procedura (dal login node):
#   module load R/4.4.2_10gcc
#   R
#   source("submit/install_packages.R")
#
# La prima volta R chiedera' di usare una libreria personale: rispondi YES
# Scegliere repository italiano quando richiesto (es. Milano o Padova)
# =============================================================================

# Repository italiano (evita di dover scegliere manualmente)
options(repos = c(CRAN = "https://cran.stat.unipd.it"))  # Padova

cat("=== PGL pipeline - R package installation ===\n")
cat("Repository:", getOption("repos"), "\n\n")

# Helper: installa solo se non presente
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing: %s ...\n", pkg))
    if (bioc) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  } else {
    cat(sprintf("OK (already installed): %s\n", pkg))
  }
}

# --- CRAN ---
cran_packages <- c(
  "BiocManager",
  "dplyr", "readr", "tidyr", "tibble",
  "ggplot2", "ggrepel", "patchwork",
  "pheatmap", "RColorBrewer",
  "openxlsx", "readxl",
  "rstatix", "ggpubr",
  "ggvenn",
  "remotes"
)

cat("--- CRAN packages ---\n")
for (pkg in cran_packages) install_if_missing(pkg)

# --- Bioconductor ---
bioc_packages <- c(
  "Seurat",
  "SingleR", "celldex",
  "DESeq2", "tximport",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "AnnotationHub",
  "Biobase",
  "SingleCellExperiment",  # dipendenza MuSiC
  "TOAST",                 # dipendenza MuSiC
  "ComplexHeatmap", "circlize"
)

cat("\n--- Bioconductor packages ---\n")
for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)

# --- MuSiC (GitHub) ---
cat("\n--- MuSiC (GitHub) ---\n")
install_if_missing("MuSiC")
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  remotes::install_github("xuranw/MuSiC")
}

cat("\n=== Installation complete ===\n")
cat("Verify with: sessionInfo()\n")
