#!/bin/bash

# This Bash script creates chromosome-specific genotype dosage matrices for 351 GTEx Heart Left Ventricle samples.
# Extracts selected samples and selected SNPs from the full GTEx VCF.
# Keeps only biallelic SNPs.
# Removes strand-ambiguous SNPs: A/T, T/A, C/G, and G/C.
# Removes variants with any missing genotype among the 351 samples.
# Changes variant IDs to CHR:POS:REF:ALT.
# Exports an additive genotype dosage matrix in PLINK .raw format.
# Deletes intermediate files unless KEEP_INTERMEDIATE=1.

# Input: 
#       GTEx WGS VCF: GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz
#       GTEx WGS VCF idx: GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz.tbi or .csi
#       SNP list: GTEx_EA_all_500kb_1_r2_0.8.prune.in
#       Sample list: heart_lv_vcf_sample_ids_not_in_original_EA110_and_in_vcf_351.txt

# Output: 
#       chr<1-22>_HLV351_dosage_nomiss.raw

set -euo pipefail

ROOT="/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS"
PLINK2="$ROOT/plink2"
VCF="/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/PopulationGenomics/GTEx/WGS_V10/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
SNPLIST="$ROOT/IntermediateResults/WGS/WGS_pruned_by_chr_500kb_1_r2_0.8/GTEx_EA_all_500kb_1_r2_0.8.prune.in"
SAMPLE_LIST="$ROOT/IntermediateResults/heart_lv_vcf_sample_ids_not_in_original_EA110_and_in_vcf_351.txt"


OUTDIR="$ROOT/IntermediateResults/WGS/by_chr_heart_left_ventricle_351_samples_variant_id"
TMP_PREFIX="$OUTDIR/HLV351_selected_snps_tmp"


mkdir -p "$OUTDIR"
"$PLINK2" \
  --vcf "$VCF" \
  --keep "$SAMPLE_LIST" \
  --extract "$SNPLIST" \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0 \
  --set-all-var-ids '@:#:$r:$a' \
  --make-pgen \
  --threads 8 \
  --out "$TMP_PREFIX"


# Export one dosage matrix per chromosome
for CHR in {1..22}; do

  echo "Exporting chromosome ${CHR}"

  "$PLINK2" \
    --pfile "$TMP_PREFIX" \
    --chr "$CHR" \
    --export A \
    --threads 8 \
    --out "$OUTDIR/chr${CHR}_HLV351_dosage_nomiss_variant_id"

done

# Delete the temporary combined PLINK dataset
rm -f \
  "${TMP_PREFIX}.pgen" \
  "${TMP_PREFIX}.pvar" \
  "${TMP_PREFIX}.psam" \
  "${TMP_PREFIX}.log"

echo "Finished exporting chromosomes 1-22."