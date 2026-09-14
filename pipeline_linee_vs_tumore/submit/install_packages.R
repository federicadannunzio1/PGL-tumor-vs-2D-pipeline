# =============================================================================
# install_packages.R
# Installa tutti i pacchetti R necessari nella home del cluster
# Esegui UNA SOLA VOLTA prima di lanciare la pipeline:
#   module load R/4.4.2_10gcc
#   Rscript submit/install_packages.R
# =============================================================================

# Directory di installazione locale (non richiede permessi di root)
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- "~/R/library"
dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)
.libPaths(lib_path)

cat("Installing to:", lib_path, "\n\n")

# CRAN
cran_packages <- c(
  "BiocManager",
  "dplyr", "readr", "tidyr", "tibble",
  "ggplot2", "ggrepel", "patchwork",
  "pheatmap", "RColorBrewer",
  "openxlsx", "readxl",
  "rstatix", "ggpubr",
  "ggvenn"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, lib = lib_path,
                     repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  } else {
    cat("OK (already installed):", pkg, "\n")
  }
}

# Bioconductor
bioc_packages <- c(
  "Seurat",
  "SingleR", "celldex",
  "DESeq2", "tximport",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "AnnotationHub",
  "Biobase",
  "ComplexHeatmap", "circlize"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "(Bioconductor)...\n")
    BiocManager::install(pkg, lib = lib_path,
                         ask = FALSE, update = FALSE)
  } else {
    cat("OK (already installed):", pkg, "\n")
  }
}

# MuSiC (GitHub - non su CRAN/Bioconductor)
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  cat("Installing MuSiC from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = lib_path,
                     repos = "https://cloud.r-project.org")
  }
  remotes::install_github("xuranw/MuSiC", lib = lib_path)
} else {
  cat("OK (already installed): MuSiC\n")
}

cat("\n=== Installation complete ===\n")
cat("Loaded library path:", lib_path, "\n")
