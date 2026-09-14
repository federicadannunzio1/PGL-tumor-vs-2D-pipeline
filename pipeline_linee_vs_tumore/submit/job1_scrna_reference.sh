#!/bin/bash
# =============================================================================
# job1_scrna_reference.sh
# Script 01: caricamento Seurat integrato (7.9 GB) + reference per MuSiC
# Partizione bigmem - nodo ad alta memoria
# =============================================================================

#SBATCH --job-name=PGL_01_scrna
#SBATCH --output=/lustre/home/gfiscon/projects/PGL/logs/job1_scrna_%j.out
#SBATCH --error=/lustre/home/gfiscon/projects/PGL/logs/job1_scrna_%j.err
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federica.dannunzio@uniroma1.it

echo "================================================"
echo "Job:   $SLURM_JOB_ID - 01_scrna_reference.R"
echo "Node:  $SLURMD_NODENAME"
echo "Start: $(date)"
echo "Mem:   120G | CPUs: $SLURM_CPUS_PER_TASK"
echo "================================================"

module load R/4.4.2_10gcc

cd /lustre/home/gfiscon/projects/PGL/pipeline || exit 1

echo "Running 01_scrna_reference.R..."
Rscript 01_scrna_reference.R

EXIT_CODE=$?
echo "Finished: $(date) | Exit code: $EXIT_CODE"
exit $EXIT_CODE
