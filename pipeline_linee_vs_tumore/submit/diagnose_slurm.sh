#!/bin/bash
# =============================================================================
# diagnose_slurm.sh
# Script di diagnostica per capire perche i job SLURM falliscono
# Esegui: bash submit/diagnose_slurm.sh
# =============================================================================

LOGDIR="/lustre/home/gfiscon/projects/PGL/logs"
mkdir -p "$LOGDIR"

echo ""
echo "============================================="
echo "DIAGNOSTICA SLURM - $(date)"
echo "============================================="

# --- 1. Partizioni disponibili ---
echo ""
echo "=== 1. PARTIZIONI DISPONIBILI ==="
sinfo -o "%.15P %.5a %.10l %.6D %.6t %.8m" 2>/dev/null || sinfo

# --- 2. Quota ---
echo ""
echo "=== 2. QUOTA ==="
lfs quota -u gfiscon /lustre/home/ 2>/dev/null

# --- 3. Test: dss + mem=1G ---
echo ""
echo "=== 3. TEST: dss mem=1G ==="
JID1=$(sbatch --parsable --partition=dss --mem=1G --time=00:02:00 \
  --output="$LOGDIR/diag_1g_%j.log" \
  --wrap="echo TEST_1G_OK && hostname && date")
echo "Submitted job: $JID1"

# --- 4. Test: dss + mem=120G ---
echo ""
echo "=== 4. TEST: dss mem=120G ==="
JID2=$(sbatch --parsable --partition=dss --mem=120G --cpus-per-task=8 --time=00:02:00 \
  --output="$LOGDIR/diag_120g_%j.log" \
  --wrap="echo TEST_120G_OK && hostname && free -h && date")
echo "Submitted job: $JID2"

# --- 5. Test: dss + mem=120G + mail ---
echo ""
echo "=== 5. TEST: dss mem=120G + mail ==="
JID3=$(sbatch --parsable --partition=dss --mem=120G --cpus-per-task=8 --time=00:02:00 \
  --mail-type=BEGIN,END,FAIL \
  --mail-user=federica.dannunzio@uniroma1.it \
  --output="$LOGDIR/diag_120g_mail_%j.log" \
  --wrap="echo TEST_MAIL_OK && hostname && date")
echo "Submitted job: $JID3"

# --- Attesa ---
echo ""
echo "Aspetto 30 secondi che i job girino..."
sleep 30

# --- Risultati ---
echo ""
echo "=== STATO JOB ==="
for JID in $JID1 $JID2 $JID3; do
  sacct -j "$JID" --format=JobID,State,ExitCode,Reason,NodeList,Start,End -X --noheader 2>/dev/null
done

echo ""
echo "=== OUTPUT LOG ==="
echo "--- Test 1G ---"
cat "$LOGDIR/diag_1g_${JID1}.log" 2>/dev/null || echo "NESSUN LOG per job $JID1"

echo ""
echo "--- Test 120G ---"
cat "$LOGDIR/diag_120g_${JID2}.log" 2>/dev/null || echo "NESSUN LOG per job $JID2"

echo ""
echo "--- Test 120G + mail ---"
cat "$LOGDIR/diag_120g_mail_${JID3}.log" 2>/dev/null || echo "NESSUN LOG per job $JID3"

echo ""
echo "============================================="
echo "DIAGNOSTICA COMPLETATA"
echo "============================================="
