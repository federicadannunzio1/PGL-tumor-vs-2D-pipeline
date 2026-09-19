#!/bin/bash
# =============================================================================
# job3_broader_deconv.sh
# Script 06-07: annotazione ampia + deconvoluzione broader (no TAM)
#
# Submit from pipeline_linee_vs_tumore/:
#   sbatch submit/job3_broader_deconv.sh
# =============================================================================

#SBATCH --job-name=PGL_06_07
#SBATCH --output=/lustre/home/gfiscon/projects/PGL/logs/job3_broader_%j.log
#SBATCH --error=/lustre/home/gfiscon/projects/PGL/logs/job3_broader_%j.log
#SBATCH --partition=dss
#SBATCH --exclude=cn2d,cn2c
#SBATCH --time=10:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federica.dannunzio@uniroma1.it

# --------------------------------------------------------------------------
# Conda
# --------------------------------------------------------------------------
source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
LOG_DIR="/lustre/home/gfiscon/projects/PGL/logs"
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------
echo "============================================="
echo "PGL pipeline — Job 06-07 (broader) — SLURM $SLURM_JOB_ID"
echo "  Date:  $(date)"
echo "  Node:  $SLURMD_NODENAME"
echo "  CPUs:  $SLURM_CPUS_PER_TASK"
echo "  RAM:   ${SLURM_MEM_PER_NODE}MB"
echo "  Dir:   $PROJECT_DIR"
echo "============================================="
R --version | head -1

# --------------------------------------------------------------------------
# Threading
# --------------------------------------------------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --------------------------------------------------------------------------
# Run scripts sequentially
# --------------------------------------------------------------------------
cd "$PROJECT_DIR" || exit 1

run_script() {
  local script=$1
  echo ""
  echo "--- Running: $script [$(date)] ---"
  Rscript "$script"
  local code=$?
  if [ $code -ne 0 ]; then
    echo "ERRORE: $script fallito (exit $code)"
    exit $code
  fi
  echo "--- Done: $script ---"
}

run_script "06_broader_annotation.R"
run_script "07_deconvolution_broader.R"

echo ""
echo "============================================="
echo "Job 06-07 completato con successo."
echo "Results: /lustre/home/gfiscon/projects/PGL/results/deconvolution_broader/"
echo "Date: $(date)"
echo "============================================="

exit 0
