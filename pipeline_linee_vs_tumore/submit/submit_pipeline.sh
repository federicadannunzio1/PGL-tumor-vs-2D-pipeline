#!/bin/bash
# =============================================================================
# submit_pipeline.sh
# Lancia la pipeline completa in due job SLURM con dipendenza:
#   Job 1 (bigmem): script 01 - carica Seurat 7.9 GB
#   Job 2 (dss):    script 02-05 - parte solo dopo che Job 1 e' completato
#
# Usage:
#   cd /lustre/home/gfiscon/projects/PGL/pipeline
#   bash submit/submit_pipeline.sh
# =============================================================================

PIPELINE_DIR="/lustre/home/gfiscon/projects/PGL/pipeline"
LOG_DIR="/lustre/home/gfiscon/projects/PGL/logs"
mkdir -p "$LOG_DIR"

echo "Submitting PGL pipeline..."

# Job 1: script 01 su bigmem (Seurat)
JOB1=$(sbatch --parsable \
  "${PIPELINE_DIR}/submit/job1_scrna_reference.sh")

echo "  Job 1 submitted: $JOB1 (01_scrna_reference.R on bigmem)"

# Job 2: script 02-05 su dss, parte solo se Job 1 e' OK
JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  "${PIPELINE_DIR}/submit/job2_bulk_degs_deconv.sh")

echo "  Job 2 submitted: $JOB2 (02-05 on dss, depends on job $JOB1)"
echo ""
echo "Monitor with:"
echo "  squeue -u gfiscon"
echo "  tail -f ${LOG_DIR}/job1_scrna_${JOB1}.out"
