# =============================================================================
# check_packages.R
# Verifica che tutti i pacchetti necessari alla pipeline siano installati
# Eseguire sul cluster prima di lanciare i job SLURM:
#   Rscript check_packages.R
# =============================================================================

required <- list(
  # 01_scrna_reference.R
  "Seurat"          = "Seurat",
  "MuSiC"           = "MuSiC",
  "Biobase"         = "Biobase",
  "patchwork"       = "patchwork",

  # 02_bulk_preprocessing.R
  "tximport"        = "tximport",
  "DESeq2"          = "DESeq2",
  "pheatmap"        = "pheatmap",
  "RColorBrewer"    = "RColorBrewer",

  # 03_degs.R
  "ggrepel"         = "ggrepel",
  "org.Hs.eg.db"    = "org.Hs.eg.db",
  "AnnotationDbi"   = "AnnotationDbi",

  # 04_deconvolution.R
  "tidyr"           = "tidyr",

  # 05_integration_figures.R
  "ggvenn"          = "ggvenn",
  "ComplexHeatmap"  = "ComplexHeatmap",
  "circlize"        = "circlize",

  # Usati ovunque
  "dplyr"           = "dplyr",
  "readr"           = "readr",
  "ggplot2"         = "ggplot2"
)

optional <- list(
  "clusterProfiler" = "clusterProfiler"   # GSEA — saltato se mancante
)

cat("=== Verifica pacchetti richiesti ===\n\n")
missing_req <- c()
for (pkg in names(required)) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-20s %s\n", pkg, ifelse(ok, "[OK]", "[MANCANTE]")))
  if (!ok) missing_req <- c(missing_req, pkg)
}

cat("\n=== Pacchetti opzionali ===\n\n")
for (pkg in names(optional)) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-20s %s\n", pkg,
              ifelse(ok, "[OK]", "[non installato - sezione GSEA verra' saltata]")))
}

cat("\n")
if (length(missing_req) == 0) {
  cat("OK: tutti i pacchetti richiesti sono disponibili. Puoi lanciare i job.\n")
} else {
  cat(sprintf("ATTENZIONE: %d pacchett%s mancant%s:\n",
              length(missing_req),
              ifelse(length(missing_req) == 1, "o", "i"),
              ifelse(length(missing_req) == 1, "e", "i")))
  for (p in missing_req) cat(sprintf("  - %s\n", p))
  cat("\nInstalla i mancanti prima di lanciare la pipeline.\n")
  quit(status = 1)
}
