# Download phased haplotypes (0/1) per chromosome using terminal


setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')

library(data.table)
simulation_dir <- "Data/Simulation"

vcf_file <- file.path(simulation_dir, "chr22_20_21Mb.vcf.gz")
out_file <- file.path(simulation_dir, "chr22_20_21Mb_1million_pseudo_subjects.csv")


# narrow to a small region


system(paste(
  "bcftools view",
  "-r 22:20000000-21000000",
  "-m2 -M2",
  "-v snps",
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
  "-Oz",
  "-o", shQuote(vcf_file)
))

system(paste( "bcftools index", shQuote(vcf_file)))

# create 100,000 pseudo subjects

set.seed(123)

n_pseudo   <- 100000L
batch_size <- 10000L

# ------------------------------------------------------------
# Read variant IDs and phased genotypes
# ------------------------------------------------------------

vcf <- data.table::fread(cmd = paste("bcftools query -f '%CHROM\\t%POS\\t%ID[\\t%GT]\\n'",
                                     shQuote(path.expand(vcf_file))), header = FALSE)
variant_info <- vcf[, 1:3]
setnames(variant_info, c("chr", "pos", "rsid"))

GT <- as.matrix(vcf[, -(1:3)])


# ------------------------------------------------------------
# Split phased genotypes into individual haplotypes
#
# Example:
#   0|1 -> haplotype 1 = 0, haplotype 2 = 1
# ------------------------------------------------------------

hap1 <- matrix( as.integer(substr(GT, 1, 1)), nrow = nrow(GT), ncol = ncol(GT))
hap2 <- matrix(as.integer(substr(GT, 3, 3)), nrow = nrow(GT), ncol = ncol(GT))


# ------------------------------------------------------------
# Haplotype pool
#
# rows    = haplotypes
# columns = SNPs
# ------------------------------------------------------------

X_pool <- rbind( t(hap1), t(hap2))


# Use rsID as SNP name; if rsID is missing, use chr:position
snp_names <- variant_info$rsid

missing_id <- is.na(snp_names) | snp_names == "."

snp_names[missing_id] <- paste0(variant_info$chr[missing_id], ":", variant_info$pos[missing_id])

colnames(X_pool) <- make.unique(snp_names)


# ------------------------------------------------------------
# Check haplotype pool
# ------------------------------------------------------------

dim(X_pool)
table(X_pool, useNA = "ifany")

# Keep SNPs with MAF >= 0.005
af  <- colMeans(X_pool, na.rm = TRUE)
maf <- pmin(af, 1 - af)
X_pool <- X_pool[, maf >= 0.025, drop = FALSE]

dim(X_pool)
# [1] 5008 2567
# ------------------------------------------------------------
# Generate 1,000,000 pseudo-subjects in batches
#
# Each subject receives two randomly sampled haplotypes:
#
# genotype = haplotype_1 + haplotype_2
#
# giving genotype dosage 0 / 1 / 2
# ------------------------------------------------------------

if (file.exists(out_file)) { file.remove(out_file)}

n_batches <- ceiling(n_pseudo / batch_size)

for (b in seq_len(n_batches)) {

  start_i <- (b - 1L) * batch_size + 1L
  end_i   <- min(b * batch_size, n_pseudo)

  n_this <- end_i - start_i + 1L

  message(
    "Batch ", b, "/", n_batches,
    ": subjects ", start_i, "-", end_i
  )


  # Randomly select two haplotypes for every pseudo-subject
  h1_idx <- sample(
    seq_len(nrow(X_pool)),
    n_this,
    replace = TRUE
  )

  h2_idx <- sample(
    seq_len(nrow(X_pool)),
    n_this,
    replace = TRUE
  )


  # Convert paired haplotypes to diploid 0/1/2 genotypes
  X <- X_pool[h1_idx, , drop = FALSE] +
    X_pool[h2_idx, , drop = FALSE]


  # Add pseudo-subject ID
  X_out <- as.data.table(X)

  X_out[
    ,
    IID := sprintf("pseudo_%07d", start_i:end_i)
  ]

  setcolorder(
    X_out,
    c("IID", setdiff(names(X_out), "IID"))
  )


  # Append each batch to the same CSV
  fwrite(
    X_out,
    file = out_file,
    append = b > 1L,
    col.names = b == 1L
  )

  rm(X, X_out)
  gc()
}


message("Finished.")
message("Saved to: ", out_file)
