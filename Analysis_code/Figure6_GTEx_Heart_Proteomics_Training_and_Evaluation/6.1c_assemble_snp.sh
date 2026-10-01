#!/bin/bash
set -euo pipefail

PAPER_DIR="${PAPER_ctOWAS_DIR:-/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS}"

PRUNE_KB="${PWAS_PRUNE_KB:-500}"
PRUNE_STEP="${PWAS_PRUNE_STEP:-1}"
PRUNE_R2="${PWAS_PRUNE_R2:-0.8}"
TAG="${PRUNE_KB}kb_${PRUNE_STEP}_r2_${PRUNE_R2}"

OUTDIR="${PWAS_MODERATE_PRUNED_DIR:-$PAPER_DIR/IntermediateResults/WGS/WGS_pruned_by_chr_${TAG}}"
OUTFILE="$OUTDIR/GTEx_EA_all_${TAG}.prune.in"

check_nonempty() {
  if [ ! -s "$1" ]; then
    echo "ERROR: missing or empty file: $1"
    exit 1
  fi
}

mkdir -p "$OUTDIR"

echo "Checking chromosome-specific prune.in files..."

for CHR in {1..22}; do
  check_nonempty "$OUTDIR/GTEx_EA_chr${CHR}_${TAG}.prune.in"
done

echo "Combining all prune.in files..."

cat "$OUTDIR"/GTEx_EA_chr{1..22}_"${TAG}".prune.in > "$OUTFILE"

check_nonempty "$OUTFILE"

echo "Combined prune list saved to:"
echo "$OUTFILE"