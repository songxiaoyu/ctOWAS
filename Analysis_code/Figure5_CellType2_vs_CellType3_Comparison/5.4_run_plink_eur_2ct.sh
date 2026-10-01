#!/bin/bash

set -euo pipefail

# PLINK2 path
PLINK2='/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS/plink2'


# Base paths
PAPER_DIR="${PAPER_ctOWAS_DIR:-/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS}"
DATA_DIR="${PAPER_DIR}/Data/plink_snplist_by_gene"
GWAS_ID_DIR="${PAPER_DIR}/Results/bcac2020_filtered_id"

mkdir -p "${GWAS_ID_DIR}"

for chr in {1..22}; do
  echo "Processing chr${chr}..."

  # Decompress files if needed (only once per chr)
  # Added -f to force overwrite if file exists to prevent errors on re-runs
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
      "${PLINK2}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
      "${PLINK2}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  # Extract selected SNPs and directly export additive dosages
  "${PLINK2}" \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/bcac2020_filtered_chr${chr}_gwas_id_pi2_comparison.txt" \
    --export A \
    --write-snplist \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_pi2_comparison"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."
