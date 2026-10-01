#!/bin/bash

# Build 1000Genome EUR LD reference files for the HERMES PWAS SNP lists.
# Run after 2_HERMES_prepare_data_pwas.R.

set -euo pipefail

PAPER_DIR="${PAPER_ctOWAS_DIR:-/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS}"
# 1000 Genomes phase 3 (hg38) reference panel, per-chromosome pgen/pvar/psam.
# Lives entirely inside Paper_PWAS -- no other project's folder is needed.
DATA_DIR="${PAPER_DIR}/Data/plink_snplist_by_gene"
# Which trained model this run is for: "ctOWAS" or "MiXcan" (must match the
# HERMES_MODEL_TAG used in 7.2_HERMES_prepare_data_pwas.R).
MODEL_TAG="${HERMES_MODEL_TAG:-ctOWAS}"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_pwas/hermes_workspace_${MODEL_TAG}"
GWAS_ID_DIR="${WORKSPACE_DIR}/hermes_filtered_id"
EUR_SAMPLES="${PAPER_DIR}/Data/1000Genome/eur_ids.txt"
CHR_LIST="$(echo {1..22})"

PLINK="${PLINK:-${PAPER_DIR}/plink2}"

for chr in ${CHR_LIST}; do
  echo "Processing chr${chr}..."

  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
    "${PLINK}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
    "${PLINK}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  "${PLINK}" \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/hermes_filtered_chr${chr}_gwas_id_pwas.txt" \
    --keep "${EUR_SAMPLES}" \
    --make-bed \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_hermes_pwas"

  "${PLINK}" \
    --bfile "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_hermes_pwas" \
    --export A \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_012_hermes_pwas"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."
