#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------
REMOTE_USER="nmille39"
REMOTE_HOST="medicinebow.arcc.uwyo.edu"

REMOTE_ROOT="/project/galacticbulge/MESA/MESA_Colors_Tests"
LOCAL_ROOT="$HOME/MESA/MESA_Colors_Tests"

# -----------------------------
# Sanity checks
# -----------------------------
if ! ssh -o BatchMode=yes -o ConnectTimeout=5 \
    "${REMOTE_USER}@${REMOTE_HOST}" 'exit 0' >/dev/null 2>&1; then
    echo "ERROR: cannot reach Medicine Bow"
    exit 1
fi

mkdir -p "$LOCAL_ROOT"

# -----------------------------
# Rsync rules
# -----------------------------
# Pull only data products:
#   - include all directories
#   - include LOGS directories and contents
#   - exclude everything else
#
# If you also want SED output, uncomment the SED include line.
# -----------------------------
rsync -avz \
  -e ssh \
  --prune-empty-dirs \
  --include '*/' \
  --include '*/LOGS*/***' \
  --include '*/LOGS_*/***' \
  --include '*/SED/***' \
  --exclude '*' \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/" \
  "${LOCAL_ROOT}/"

echo "Data sync complete."
