# =============================================================================
# 09_low_expression_analysis.R
# Analisi di geni attesi a bassa espressione sia nel tumore che nelle linee 2D
# Geni target: KREMEN2, HIF1A, TH, PHGDH, CDH1
#
# Input:  scRNA integrato (Seurat), TPM matrix (bulk), DEG table, metadata
# Output: violin/dot/feature plot (scRNA), boxplot bulk, summary CSV
# =============================================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Risolvi conflitto select: AnnotationDbi maschera dplyr::select
select <- dplyr::select

# -----------------------------------------------------------------------------
# CHECK INPUT
# -----------------------------------------------------------------------------
check_inputs(
  SCRNA_INTEGRATED_RDS,
  file.path(RESULTS_BULK, "tpm_matrix.csv"),
  SAMPLE_METADATA
)

# -----------------------------------------------------------------------------
# OUTPUT DIRECTORY
# -----------------------------------------------------------------------------
RESULTS_LOWEXPR <- file.path(RESULTS_DIR, "low_expression")
dir.create(RESULTS_LOWEXPR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# GENI TARGET
# -----------------------------------------------------------------------------
target_genes <- c("KREMEN2", "HIF1A", "TH", "PHGDH", "CDH1")

# =============================================================================
# PARTE A — Analisi scRNA-seq
# =============================================================================
message("\n--- PARTE A: scRNA-seq ---")

# Caricamento Seurat
message("Caricamento oggetto Seurat...")
seu <- readRDS(SCRNA_INTEGRATED_RDS)

seu$celltype_pasquale <- as.character(Idents(seu))
Idents(seu) <- "celltype_pasquale"

# Check geni disponibili
available_genes <- target_genes[target_genes %in% rownames(seu)]
missing_genes   <- target_genes[!target_genes %in% rownames(seu)]

if (length(missing_genes) > 0) {
  message("Geni non trovati nella matrice scRNA: ", paste(missing_genes, collapse = ", "))
}
message("Geni disponibili per scRNA: ", paste(available_genes, collapse = ", "))

if (length(available_genes) > 0) {

  # --- VlnPlot per ogni gene ---
  message("Creazione VlnPlot...")
  vln_list <- lapply(available_genes, function(g) {
    VlnPlot(seu, features = g, pt.size = 0) +
      THEME_PGL +
      ggtitle(g)
  })

  pdf(file.path(RESULTS_LOWEXPR, "vlnplot_per_gene.pdf"),
      width = 12, height = 4 * length(available_genes))
  print(wrap_plots(vln_list, ncol = 1))
  dev.off()
  message("Salvato: vlnplot_per_gene.pdf")

  # --- DotPlot ---
  message("Creazione DotPlot...")
  pdf(file.path(RESULTS_LOWEXPR, "dotplot_target_genes.pdf"),
      width = 10, height = 6)
  print(
    DotPlot(seu, features = available_genes) +
      RotatedAxis() +
      THEME_PGL +
      ggtitle("Target gene expression by cell type")
  )
  dev.off()
  message("Salvato: dotplot_target_genes.pdf")

  # --- FeaturePlot UMAP per ogni gene ---
  message("Creazione FeaturePlot...")
  feat_list <- lapply(available_genes, function(g) {
    FeaturePlot(seu, features = g, reduction = "umap") +
      THEME_PGL +
      ggtitle(g)
  })

  pdf(file.path(RESULTS_LOWEXPR, "featureplot_umap.pdf"),
      width = 10, height = 4 * ceiling(length(available_genes) / 2))
  print(wrap_plots(feat_list, ncol = 2))
  dev.off()
  message("Salvato: featureplot_umap.pdf")

} else {
  message("Nessun gene target trovato nella matrice scRNA. Salto Parte A.")
}

# Libera memoria
rm(seu); gc()

# =============================================================================
# PARTE B — Confronto bulk RNA-seq
# =============================================================================
message("\n--- PARTE B: Bulk RNA-seq ---")

# Caricamento TPM
tpm_mat <- as.data.frame(read_csv(file.path(RESULTS_BULK, "tpm_matrix.csv"),
                                  show_col_types = FALSE))
rownames(tpm_mat) <- tpm_mat[[1]]
tpm_mat <- tpm_mat[, -1]

# Strip versione ENSEMBL e rimuovi duplicati
stripped <- sub("\\..*", "", rownames(tpm_mat))
tpm_mat <- tpm_mat[!duplicated(stripped), ]
rownames(tpm_mat) <- stripped[!duplicated(stripped)]

# Caricamento metadata
meta <- read_csv(SAMPLE_METADATA, show_col_types = FALSE)

# Caricamento annotazione geni (ENSEMBL -> symbol) dalla tabella DEG
degs_full <- read_csv(file.path(RESULTS_DEGS, "degs_full_table.csv"),
                      show_col_types = FALSE)

# Costruisci mapping ENSEMBL -> symbol
gene_map <- degs_full %>%
  select(ensembl_id = gene_id, symbol = gene_symbol) %>%
  distinct()
gene_map$ensembl_id <- sub("\\..*", "", gene_map$ensembl_id)

# Trova ENSEMBL ID per ogni gene target
target_map <- gene_map %>%
  filter(symbol %in% target_genes)

message("Mapping geni target:")
print(target_map)

if (nrow(target_map) > 0) {

  # Estrai TPM per geni target
  target_ensembl <- target_map$ensembl_id
  found_in_tpm   <- target_ensembl[target_ensembl %in% rownames(tpm_mat)]

  if (length(found_in_tpm) > 0) {

    tpm_target <- tpm_mat[found_in_tpm, , drop = FALSE]
    tpm_target$ensembl_id <- rownames(tpm_target)

    tpm_long <- tpm_target %>%
      pivot_longer(-ensembl_id, names_to = "sample_id", values_to = "tpm") %>%
      left_join(target_map, by = "ensembl_id") %>%
      left_join(meta %>% select(sample_id, condition), by = "sample_id") %>%
      mutate(log2_tpm = log2(tpm + 1))

    # Boxplot
    p_bulk <- ggplot(tpm_long, aes(x = condition, y = log2_tpm, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
      facet_wrap(~ symbol, scales = "free_y") +
      scale_fill_manual(values = c("tumor" = "#E64B35", "2D" = "#4DBBD5")) +
      labs(x = "Condition", y = "log2(TPM + 1)",
           title = "Target gene expression \u2014 Bulk RNA-seq") +
      THEME_PGL +
      theme(legend.position = "none")

    pdf(file.path(RESULTS_LOWEXPR, "boxplot_bulk_target_genes.pdf"),
        width = 10, height = 8)
    print(p_bulk)
    dev.off()
    message("Salvato: boxplot_bulk_target_genes.pdf")

  } else {
    message("Nessun gene target trovato nella matrice TPM.")
  }
} else {
  message("Nessun mapping ENSEMBL trovato per i geni target.")
}

# =============================================================================
# PARTE C — Tabella riassuntiva
# =============================================================================
message("\n--- PARTE C: Summary table ---")

if (nrow(target_map) > 0 && length(found_in_tpm) > 0) {

  # Media TPM per condizione
  summary_tpm <- tpm_long %>%
    group_by(symbol, condition) %>%
    summarise(mean_tpm = mean(tpm, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = mean_tpm,
                names_prefix = "mean_tpm_")

  # LFC e padj da DESeq2
  degs_summary <- degs_full %>%
    filter(gene_symbol %in% target_genes) %>%
    select(symbol = gene_symbol, log2FC = log2FoldChange, padj) %>%
    distinct(symbol, .keep_all = TRUE)

  summary_table <- summary_tpm %>%
    left_join(degs_summary, by = "symbol")

  write_csv(summary_table,
            file.path(RESULTS_LOWEXPR, "summary_low_expression_genes.csv"))
  message("Salvato: summary_low_expression_genes.csv")

  message("\n--- Tabella riassuntiva ---")
  print(as.data.frame(summary_table))

} else {
  message("Impossibile creare summary: mapping geni non disponibile.")
}

message("\n=== 09_low_expression_analysis.R completato ===")
