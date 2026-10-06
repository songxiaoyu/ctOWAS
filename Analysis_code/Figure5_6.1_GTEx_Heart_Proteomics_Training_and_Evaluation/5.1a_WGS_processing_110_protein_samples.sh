#!/bin/bash

# This Bash script creates chromosome-specific genotype dosage matrices for 351 GTEx Heart Left Ventricle samples.
# Extracts selected samples and selected SNPs from the full GTEx VCF.
# Keeps only biallelic SNPs.
# Removes strand-ambiguous SNPs: A/T, T/A, C/G, and G/C.
# Removes variants with any missing genotype.
# Changes variant IDs to CHR:POS:REF:ALT.
# Exports an additive genotype dosage matrix in PLINK .raw format.
# Deletes intermediate files unless KEEP_INTERMEDIATE=1.

# Input: 
#       GTEx WGS VCF: GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz
#       GTEx WGS VCF idx: GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz.tbi or .csi
#       SNP list: 1000g_b38_snpIDs.txt
#       Sample list: GTEx_EA_samples.txt

# Output: 
#       chr<1-22>_HLV351_dosage_nomiss.raw


set -euo pipefail

PLINK="./plink2"
# /Users/zhusinan/Library/CloudStorage/Dropbox/Paper_ctOWAS/New generated files/codes
VCF="$HOME/Documents/vcf/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
SAMPLES_PLINK="$HOME/Documents/vcf/GTEx_EA_samples.txt"
SAMPLES_BCF="$HOME/Documents/vcf/GTEx_EA_samples_bcftools.txt"
SNPLIST="$HOME/Documents/vcf/1000g_b38_snpIDs.txt"
OUTDIR="$HOME/Documents/vcf/by_chr_nomiss"

mkdir -p "$OUTDIR"

echo "Preparing bcftools sample file..."
awk '{print $1}' "$SAMPLES_PLINK" | tr -d '\r' > "$SAMPLES_BCF"

check_nonempty() {
  if [ ! -s "$1" ]; then
    echo "ERROR: missing or empty file: $1"
    exit 1
  fi
}

for CHR in {1..22}; do
  echo "===== Processing chr${CHR} ====="

  CHR_VCF="$OUTDIR/chr${CHR}_GTEx_EA.vcf.gz"
  CHR_TBI="$CHR_VCF.tbi"
  PLINK_VCF="$OUTDIR/chr${CHR}_filtered.vcf"
  FINAL_VCF="$OUTDIR/GTEx_EA_chr${CHR}_filtered.vcf.gz"
  FINAL_TBI="$FINAL_VCF.tbi"
  DOSAGE_RAW="$OUTDIR/chr${CHR}_dosage_nomiss.raw"


  # 1. Extract chromosome + EA samples
  if [ -s "$CHR_VCF" ]; then
    echo "Step 1 exists, skipping: $CHR_VCF"
  else
    echo "Step 1: extracting chr${CHR}"
    bcftools view \
      -r chr${CHR} \
      -S "$SAMPLES_BCF" \
      -I \
      -Oz \
      -o "$CHR_VCF" \
      "$VCF"
  fi
  check_nonempty "$CHR_VCF"

  # 2. Index chr VCF
  if [ -s "$CHR_TBI" ]; then
    echo "Step 2 exists, skipping: $CHR_TBI"
  else
    echo "Step 2: indexing chr${CHR} VCF"
    bcftools index -t "$CHR_VCF"
  fi
  check_nonempty "$CHR_TBI"

  # 3. Extract SNP list with PLINK
  if [ -s "$PLINK_VCF" ]; then
    echo "Step 3 exists, skipping: $PLINK_VCF"
  else
    echo "Step 3: PLINK extract SNPs chr${CHR}"
    "$PLINK" \
      --vcf "$CHR_VCF" \
      --extract "$SNPLIST" \
      --export vcf \
      --force-intersect \
      --out "$OUTDIR/chr${CHR}_filtered"
  fi
  check_nonempty "$PLINK_VCF"

  # 4. Keep only single-base SNPs, remove strand-ambiguous SNPs, and get rid of missingness
  if [ -s "$FINAL_VCF" ]; then
    echo "Step 4 exists, skipping: $FINAL_VCF"
  else
    echo "Step 4: filtering SNPs chr${CHR}"

    bcftools view \
      -i '(F_MISSING=0) && (STRLEN(REF)==1) && (STRLEN(ALT)==1) && (ALT!="A" || REF!="T") && (ALT!="T" || REF!="A") && (ALT!="C" || REF!="G") && (ALT!="G" || REF!="C")' \
      -Oz \
      -o "$FINAL_VCF" \
      "$PLINK_VCF"
  fi
  check_nonempty "$FINAL_VCF"

  # 4b. Index final filtered VCF
  if [ -s "$FINAL_TBI" ]; then
    echo "Step 4b exists, skipping: $FINAL_TBI"
  else
    echo "Step 4b: indexing final filtered VCF chr${CHR}"
    bcftools index -t "$FINAL_VCF"
  fi
  check_nonempty "$FINAL_TBI"

  # 4c. Convert to a PLINK2 pgen fileset with standardized CHR:POS:REF:ALT
  #     variant IDs. This is the "*_nomiss" fileset that 6.1b (LD pruning)
  #     expects as input, and it also makes the 110-sample training variant
  #     IDs use the same CHR:POS:REF:ALT scheme as the 351-sample prediction
  #     cohort (6.1d), so SNPs line up directly by ID instead of relying on
  #     the rsID-to-position fallback lookup in helper_function.R.
  NOMISS_PREFIX="$OUTDIR/GTEx_EA_chr${CHR}_nomiss"
  if [ -s "${NOMISS_PREFIX}.pgen" ]; then
    echo "Step 4c exists, skipping: ${NOMISS_PREFIX}.pgen"
  else
    echo "Step 4c: exporting pgen fileset chr${CHR}"
    "$PLINK" \
      --vcf "$FINAL_VCF" \
      --set-all-var-ids '@:#:$r:$a' \
      --make-pgen \
      --out "$NOMISS_PREFIX"
  fi
  check_nonempty "${NOMISS_PREFIX}.pgen"

  # 5. Export dosage/additive genotype matrix
  if [ -s "$DOSAGE_RAW" ]; then
    echo "Step 5 exists, skipping: $DOSAGE_RAW"
  else
    echo "Step 5: exporting dosage chr${CHR}"
    "$PLINK" \
      --pfile "$NOMISS_PREFIX" \
      --export A \
      --out "$OUTDIR/chr${CHR}_dosage_nomiss"
  fi
  check_nonempty "$DOSAGE_RAW"

  echo "===== Finished chr${CHR} ====="
done