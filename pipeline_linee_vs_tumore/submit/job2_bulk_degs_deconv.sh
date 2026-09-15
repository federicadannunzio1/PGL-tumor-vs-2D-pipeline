#!/bin/bash
# =============================================================================
# job2_bulk_degs_deconv.sh
# Script 02-05: preprocessing bulk, DEGs, deconvoluzione, figure
# Partizione dss - memoria standard
# Questo job parte SOLO se job1 e' completato con successo
# =============================================================================

#SBATCH --job-name=PGL_02_05
#SBATCH --output=/lustre/home/gfiscon/projects/PGL/logs/job2_pipeline_%j.out
#SBATCH --error=/lustre/home/gfiscon/projects/PGL/logs/job2_pipeline_%j.err
#SBATCH --partition=dss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federica.dannunzio@uniroma1.it

echo "================================================"
echo "Job:   $SLURM_JOB_ID - scripts 02-05"
echo "Node:  $SLURMD_NODENAME"
echo "Start: $(date)"
echo "Mem:   48G | CPUs: $SLURM_CPUS_PER_TASK"
echo "================================================"

source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

PIPELINE_DIR="/lustre/home/gfiscon/projects/PGL/pipeline/pipeline_linee_vs_tumore"
cd "$PIPELINE_DIR" || exit 1

run_script() {
  local script=$1
  echo ""
  echo "--- Running: $script [$(date)] ---"
  Rscript "$script"
  local code=$?
  if [ $code -ne 0 ]; then
    echo "ERROR: $script failed (exit $code)"
    exit $code
  fi
  echo "--- Done: $script ---"
}

run_script "02_bulk_preprocessing.R"
run_script "03_degs.R"
run_script "04_deconvolution.R"
run_script "05_integration_figures.R"

echo ""
echo "================================================"
echo "Pipeline completed: $(date)"
echo "Results: /lustre/home/gfiscon/projects/PGL/results/"
echo "================================================"
