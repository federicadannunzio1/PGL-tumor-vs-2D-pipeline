#!/bin/bash
# =============================================================================
# job1_scrna_reference.sh
# Script 01: caricamento Seurat integrato + reference per MuSiC
#
# Submit from pipeline_linee_vs_tumore/:
#   sbatch submit/job1_scrna_reference.sh
# =============================================================================

#SBATCH --job-name=PGL_01_scrna
#SBATCH --output=/lustre/home/gfiscon/projects/PGL/logs/job1_scrna_%j.log
#SBATCH --error=/lustre/home/gfiscon/projects/PGL/logs/job1_scrna_%j.log
#SBATCH --partition=dss
#SBATCH --exclude=cn2d
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
# Paths — derived from submission directory
# --------------------------------------------------------------------------
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
LOG_DIR="/lustre/home/gfiscon/projects/PGL/logs"
mkdir -p "$LOG_DIR"

# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------
echo "============================================="
echo "PGL pipeline — Job 01 — SLURM $SLURM_JOB_ID"
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
# Run
# --------------------------------------------------------------------------
cd "$PROJECT_DIR" || exit 1

echo "Running: 01_scrna_reference.R"
Rscript 01_scrna_reference.R
EXIT_CODE=$?

echo ""
echo "============================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "Job 01 completato con successo."
else
  echo "Job 01 FALLITO (exit code: $EXIT_CODE)."
  echo "Log: $LOG_DIR/job1_scrna_${SLURM_JOB_ID}.log"
fi
echo "Date: $(date)"
echo "============================================="

exit $EXIT_CODE
