#!/bin/bash
# =============================================================================
# submit_pipeline.sh
# SLURM job submission - PGL tumor vs 2D primary lines pipeline
#
# Usage:
#   cd /lustre/home/gfiscon/projects/PGL/pipeline
#   sbatch submit/submit_pipeline.sh
# =============================================================================

#SBATCH --job-name=PGL_pipeline
#SBATCH --output=/lustre/home/gfiscon/projects/PGL/logs/pipeline_%j.out
#SBATCH --error=/lustre/home/gfiscon/projects/PGL/logs/pipeline_%j.err
#SBATCH --partition=cpu           # TODO: verifica nome partizione con: sinfo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G                 # 64G per Seurat + margine per gli altri script
#SBATCH --time=12:00:00           # 12 ore (conservativo)
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federica.dannunzio@uniroma1.it

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
echo "============================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURMD_NODENAME"
echo "Start:     $(date)"
echo "============================================"

PIPELINE_DIR="/lustre/home/gfiscon/projects/PGL/pipeline"
LOG_DIR="/lustre/home/gfiscon/projects/PGL/logs"

mkdir -p "$LOG_DIR"
cd "$PIPELINE_DIR" || { echo "ERROR: cannot cd to $PIPELINE_DIR"; exit 1; }

# Carica modulo R (verifica versione disponibile con: module avail R)
module load R

echo "R version: $(R --version | head -1)"

# -----------------------------------------------------------------------------
# Funzione helper: esegui script R e controlla errori
# -----------------------------------------------------------------------------
run_script() {
  local script=$1
  echo ""
  echo "--------------------------------------------"
  echo "Running: $script  [$(date)]"
  echo "--------------------------------------------"
  Rscript "$script"
  local exit_code=$?
  if [ $exit_code -ne 0 ]; then
    echo "ERROR: $script failed with exit code $exit_code"
    exit $exit_code
  fi
  echo "Done: $script  [$(date)]"
}

# -----------------------------------------------------------------------------
# Pipeline - esecuzione sequenziale
# -----------------------------------------------------------------------------
run_script "02_bulk_preprocessing.R"
run_script "01_scrna_reference.R"    # piu' lento, carica 7.9 GB Seurat
run_script "03_degs.R"
run_script "04_deconvolution.R"
run_script "05_integration_figures.R"

# -----------------------------------------------------------------------------
# Fine
# -----------------------------------------------------------------------------
echo ""
echo "============================================"
echo "Pipeline completed: $(date)"
echo "Results: /lustre/home/gfiscon/projects/PGL/results/"
echo "============================================"
