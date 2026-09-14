#!/bin/bash
# =============================================================================
# sync_data_to_cluster.sh
# Trasferisce i dati da Google Drive al cluster Terastat
#
# PREREQUISITI:
#   1. rclone configurato con remote "gdrive" puntato al tuo Google Drive
#      (se non e' configurato: rclone config)
#   2. Oppure, per i file piccoli (SCEVAN, bulk):
#      scarica offline da Google Drive, poi usa rsync
#
# Usage (da terminale Mac):
#   bash submit/sync_data_to_cluster.sh
# =============================================================================

CLUSTER_USER="gfiscon"
CLUSTER_HOST="terastat.uniroma1.it"   # TODO: verifica hostname esatto
CLUSTER_BASE="/lustre/home/gfiscon/projects/PGL"

# Path Google Drive locale (dove sono i file se scaricati offline)
GD_BASE="/Users/federicadannunzio/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/caruana-project/PGL"

# Remote rclone (nome configurato con: rclone config)
RCLONE_REMOTE="gdrive"
RCLONE_GD_PATH="Drive condivisi/caruana-project/PGL"

echo "============================================"
echo "PGL data sync to Terastat cluster"
echo "Target: ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}"
echo "============================================"

# -----------------------------------------------------------------------------
# STEP 1: Crea struttura cartelle sul cluster
# -----------------------------------------------------------------------------
echo ""
echo "[1/5] Creating directory structure on cluster..."
ssh "${CLUSTER_USER}@${CLUSTER_HOST}" "
  mkdir -p ${CLUSTER_BASE}/pipeline/data
  mkdir -p ${CLUSTER_BASE}/pipeline/submit
  mkdir -p ${CLUSTER_BASE}/data/scrna_integrated
  mkdir -p ${CLUSTER_BASE}/data/bulk_salmon
  mkdir -p ${CLUSTER_BASE}/data/scevan/PC190
  mkdir -p ${CLUSTER_BASE}/data/scevan/PTJ173
  mkdir -p ${CLUSTER_BASE}/data/scevan/PTJ184
  mkdir -p ${CLUSTER_BASE}/data/scevan/PTJ185
  mkdir -p '${CLUSTER_BASE}/data/scevan/PV158 BIS'
  mkdir -p ${CLUSTER_BASE}/data/scevan/PV180
  mkdir -p ${CLUSTER_BASE}/data/scevan/PV181
  mkdir -p ${CLUSTER_BASE}/data/scevan/PV193
  mkdir -p ${CLUSTER_BASE}/reference
  mkdir -p ${CLUSTER_BASE}/results
  mkdir -p ${CLUSTER_BASE}/logs
  echo 'Directory structure created.'
"

# -----------------------------------------------------------------------------
# STEP 2: Trasferisci gli script R (dal git clone o direttamente)
# Opzione A: git clone sul cluster (consigliata, aggiorna con git pull)
# -----------------------------------------------------------------------------
echo ""
echo "[2/5] Syncing R scripts..."
rsync -avz --progress \
  "${GD_BASE}/analisi_fede/pipeline_linee_vs_tumore/"*.R \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/pipeline/"

rsync -avz --progress \
  "${GD_BASE}/analisi_fede/pipeline_linee_vs_tumore/submit/" \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/pipeline/submit/"

rsync -avz --progress \
  "${GD_BASE}/analisi_fede/pipeline_linee_vs_tumore/data/sample_metadata.csv" \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/pipeline/data/"

# -----------------------------------------------------------------------------
# STEP 3: Trasferisci file Salmon bulk RNA-seq (file piccoli, rsync diretto)
# NOTA: assicurati che i file siano disponibili offline in Google Drive
# -----------------------------------------------------------------------------
echo ""
echo "[3/5] Syncing Salmon bulk RNA-seq files..."
rsync -avz --progress \
  --include="*.quant.genes.sf" \
  --exclude="*.quant.sf" \
  "${GD_BASE}/analisi_pasquale/bulk RNA/RNA counts/" \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/data/bulk_salmon/"

# -----------------------------------------------------------------------------
# STEP 4: Trasferisci file SCEVAN per-sample (solo scevan_PGL.RDS)
# -----------------------------------------------------------------------------
echo ""
echo "[4/5] Syncing SCEVAN per-sample RDS files..."
for SAMPLE in PC190 PTJ173 PTJ184 PTJ185 PV180 PV181 PV193; do
  rsync -avz --progress \
    "${GD_BASE}/analisi_alessio/scevan_iterato_tutti_campioni/${SAMPLE}/scevan_PGL.RDS" \
    "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/data/scevan/${SAMPLE}/"
done

# PV158 BIS ha uno spazio nel nome
rsync -avz --progress \
  "${GD_BASE}/analisi_alessio/scevan_iterato_tutti_campioni/PV158 BIS/scevan_PGL.RDS" \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/data/scevan/PV158 BIS/"

# Marker file (se disponibile offline)
rsync -avz --progress \
  "${GD_BASE}/analisi_pasquale/scRNA PGL/results/marker_per_cluster.xlsx" \
  "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/reference/" 2>/dev/null || \
  echo "  WARNING: marker_per_cluster.xlsx non disponibile offline, skip."

# -----------------------------------------------------------------------------
# STEP 5: Seurat integrato (7.9 GB) - USA RCLONE, non rsync
# rclone va direttamente da Google Drive al cluster via SSH/sftp
# -----------------------------------------------------------------------------
echo ""
echo "[5/5] Transferring Seurat integrated object (7.9 GB) via rclone..."
echo "      This may take 20-60 minutes depending on your connection."
echo ""

# Verifica che rclone sia installato
if ! command -v rclone &> /dev/null; then
  echo "  WARNING: rclone non trovato."
  echo "  Installalo con: brew install rclone"
  echo "  Configuralo con: rclone config"
  echo ""
  echo "  In alternativa, sul cluster esegui:"
  echo "    rclone copy \"${RCLONE_REMOTE}:${RCLONE_GD_PATH}/analisi_alessio/PGL_PFE_3_integrated_x_DE_complete_annotation.Rds\" \\"
  echo "      ${CLUSTER_BASE}/data/scrna_integrated/"
else
  rclone copy \
    "${RCLONE_REMOTE}:${RCLONE_GD_PATH}/analisi_alessio/PGL_PFE_3_integrated_x_DE_complete_annotation.Rds" \
    "${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_BASE}/data/scrna_integrated/" \
    --sftp-host "${CLUSTER_HOST}" \
    --sftp-user "${CLUSTER_USER}" \
    --progress \
    --transfers 1
fi

echo ""
echo "============================================"
echo "Sync completed: $(date)"
echo "============================================"
echo ""
echo "NEXT STEPS on cluster:"
echo "  cd ${CLUSTER_BASE}/pipeline"
echo "  sbatch submit/submit_pipeline.sh"
echo "  squeue -u ${CLUSTER_USER}"
