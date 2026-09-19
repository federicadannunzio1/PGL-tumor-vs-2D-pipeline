# =============================================================================
# 08_featureplots_genes.R
# FeaturePlot e DotPlot su UMAP del Seurat integrato di Pasquale
#
# Obiettivi:
# Visualizzare l'espressione di geni di interesse (mesenchimali, WNT, S100,
# PAX, enolasi) sulle popolazioni cellulari identificate nel tumore.
#
# Input:  oggetto Seurat integrato (pasquale_s_obj_multi_finalAnnotation.rds)
# Output: PDF con FeaturePlot e DotPlot per ogni gruppo di geni
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})

check_inputs(SCRNA_INTEGRATED_RDS)

# Output directory
RESULTS_FEAT <- file.path(RESULTS_DIR, "featureplots")
dir.create(RESULTS_FEAT, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO SEURAT
# -----------------------------------------------------------------------------
message("\n--- 1. Caricamento Seurat ---")

seu <- readRDS(SCRNA_INTEGRATED_RDS)
seu$celltype_pasquale <- as.character(Idents(seu))
Idents(seu) <- "celltype_pasquale"

cat(sprintf("Cellule: %d | Geni: %d\n", ncol(seu), nrow(seu)))
cat("Tipi cellulari:\n")
print(sort(table(Idents(seu)), decreasing = TRUE))

# -----------------------------------------------------------------------------
# 2. DEFINIZIONE GRUPPI DI GENI
# -----------------------------------------------------------------------------
gene_groups <- list(
  "Mesenchymal_Stromal" = c("PDGFRA", "COL6A1", "COL6A2", "COL6A3",
                             "VIM", "GREM1"),
  "WNT_pathway"         = c("WNT5A", "WNT5B"),
  "S100_family"         = c("S100A2", "S100A10", "S100A11", "S100A13",
                             "S100A16", "S100B", "S100P", "S100A1"),
  "PAX_family"          = c("PAX8", "PAX3", "PAX8-AS1"),
  "Enolases"            = c("ENO1", "ENO2", "ENO3")
)

# -----------------------------------------------------------------------------
# 3. PER OGNI GRUPPO: FEATUREPLOT + DOTPLOT
# -----------------------------------------------------------------------------
message("\n--- 2. FeaturePlot e DotPlot per gruppo ---")

all_available_genes <- c()
gene_to_group       <- c()

for (group_name in names(gene_groups)) {
  genes     <- gene_groups[[group_name]]
  available <- genes[genes %in% rownames(seu)]
  missing   <- genes[!genes %in% rownames(seu)]

  cat(sprintf("\n--- Gruppo: %s ---\n", group_name))
  cat(sprintf("  Richiesti: %d | Disponibili: %d | Mancanti: %d\n",
              length(genes), length(available), length(missing)))
  if (length(missing) > 0) cat("  MANCANTI:", paste(missing, collapse = ", "), "\n")

  if (length(available) == 0) {
    message("  Nessun gene disponibile, salto questo gruppo.")
    next
  }

  all_available_genes <- c(all_available_genes, available)
  gene_to_group       <- c(gene_to_group, rep(group_name, length(available)))

  # FeaturePlot
  n_cols <- min(3, length(available))
  n_rows <- ceiling(length(available) / n_cols)

  fp <- FeaturePlot(
    seu, features = available,
    reduction = "umap", pt.size = 0.3, order = TRUE, ncol = n_cols
  )

  pdf(file.path(RESULTS_FEAT, sprintf("featureplot_%s.pdf", group_name)),
      width = n_cols * 5, height = n_rows * 4.5)
  print(fp)
  dev.off()

  # DotPlot
  dp <- DotPlot(seu, features = available) +
    labs(title = sprintf("DotPlot - %s", gsub("_", " ", group_name))) +
    THEME_PGL +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

  ggsave(file.path(RESULTS_FEAT, sprintf("dotplot_%s.pdf", group_name)),
         dp, width = max(6, length(available) * 0.8 + 3),
         height = 7)
}

# -----------------------------------------------------------------------------
# 4. DOTPLOT COMBINATO — TUTTI I GENI
# -----------------------------------------------------------------------------
message("\n--- 3. DotPlot combinato ---")

if (length(all_available_genes) > 0) {
  dp_all <- DotPlot(seu, features = all_available_genes) +
    labs(title = "Expression of selected genes by cell type") +
    THEME_PGL +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

  ggsave(file.path(RESULTS_FEAT, "dotplot_all_genes_combined.pdf"),
         dp_all,
         width  = max(10, length(all_available_genes) * 0.7 + 3),
         height = 8)

  write_csv(
    data.frame(gene = all_available_genes, group = gene_to_group),
    file.path(RESULTS_FEAT, "genes_plotted_summary.csv")
  )
}

message("\n=== 08_featureplots_genes.R completato ===")
message("Output in: ", RESULTS_FEAT)
