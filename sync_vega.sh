#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------

REMOTE_USER="nill"
REMOTE_HOST="192.168.0.123"   # vega
JUMP_HOST="cirrus"

REMOTE_ROOT="/home/nill/MESA/MESA_Colors_Tests"
LOCAL_ROOT="$HOME/MESA/MESA_Colors_Tests"

SSH_OPTS="-J ${JUMP_HOST}"

# -----------------------------
# Sanity checks
# -----------------------------

if ! ssh -o BatchMode=yes -o ConnectTimeout=5 $SSH_OPTS \
     "${REMOTE_USER}@${REMOTE_HOST}" 'exit 0' >/dev/null 2>&1; then
    echo "ERROR: cannot reach vega via cirrus"
    exit 1
fi

mkdir -p "$LOCAL_ROOT"

# -----------------------------
# Rsync rules
# -----------------------------
# We:
#   - include all directories
#   - include */SED/** and */LOGS/**
#   - exclude everything else
#
# This guarantees we only pull data products.
# -----------------------------

rsync -avz --delete \
  -e "ssh $SSH_OPTS" \
  --prune-empty-dirs \
  --include '*/' \
  --include '*/SED/***' \
  --include '*/LOGS/***' \
  --exclude '*' \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/" \
  "${LOCAL_ROOT}/"

echo "Data sync complete."
