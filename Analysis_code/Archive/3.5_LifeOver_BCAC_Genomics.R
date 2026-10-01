# NOTE: DRIVE Genotype is already Build 38.

rm(list=ls())
library(dplyr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

setwd('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS')
paper_dir <- getwd()
drive_geno_dir <- file.path('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_MiXcan/Data/OncoArray')
results_dir <- file.path(paper_dir, "Results", "DRIVE")

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(R.utils)

# Import the uncompressed hg19-to-hg38 chain file once
chain_file <- file.path(paper_dir, "Data", "hg19ToHg38.over.chain")
# if (!file.exists(chain_file)) R.utils::gunzip(chain_gz, destname = chain_file, overwrite = TRUE, remove = FALSE)

hg19_to_hg38_chain <- import.chain(chain_file)

liftover_summary <- vector("list", 22L)

for (chr in seq_len(22L)) {

  message("Processing chromosome ", chr, " of 22")

  input_file <- file.path(drive_geno_dir, sprintf("oncoarray_dosages_chr%02d.txt.gz", chr))
  map_file <- file.path(drive_geno_dir, sprintf("oncoarray_chr%02d_build38.txt", chr))
  unmapped_file <- file.path(drive_geno_dir, sprintf("oncoarray_chr%02d_build37_unmapped.txt", chr))

  if (!file.exists(input_file)) {
    warning("Input file unavailable: ", input_file)
    liftover_summary[[chr]] <- data.table(CHR = chr, input_variants = 0L, mapped_variants = 0L, unmapped_variants = 0L, mapping_rate = NA_real_, status = "file_missing")
    next
  }

  # Read only CHR, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT.
  # The tens of thousands of sample-dosage columns are not loaded.
  drive_variant_info <- fread(
    file = input_file,
    header = TRUE,
    select = 1:9,
    data.table = TRUE,
    check.names = FALSE
  )

  if ("#CHROM" %in% names(drive_variant_info)) setnames(drive_variant_info, "#CHROM", "CHR")
  if ("CHROM" %in% names(drive_variant_info)) setnames(drive_variant_info, "CHROM", "CHR")
  if (!"CHR" %in% names(drive_variant_info)) setnames(drive_variant_info, 1L, "CHR")

  stopifnot(all(c("CHR", "POS", "ID", "REF", "ALT") %in% names(drive_variant_info)))

  drive_variant_info[, `:=`(
    row_id = .I,
    CHR = as.character(CHR),
    POS = as.integer(POS)
  )]

  drive_gr37 <- GRanges(
    seqnames = ifelse(grepl("^chr", drive_variant_info$CHR, ignore.case = TRUE), drive_variant_info$CHR, paste0("chr", drive_variant_info$CHR)),
    ranges = IRanges(start = drive_variant_info$POS, width = 1L),
    row_id = drive_variant_info$row_id
  )

  drive_gr38_list <- liftOver(drive_gr37, hg19_to_hg38_chain)
  n_lifted <- lengths(drive_gr38_list)
  one_to_one <- which(n_lifted == 1L)

  drive_liftover_map <- data.table(
    row_id = drive_variant_info$row_id,
    n_lifted = n_lifted,
    CHR_build38 = NA_character_,
    POS_build38 = NA_integer_
  )

  if (length(one_to_one)) {
    drive_gr38 <- unlist(drive_gr38_list[one_to_one], use.names = FALSE)
    drive_liftover_map[one_to_one, `:=`(CHR_build38 = as.character(seqnames(drive_gr38)), POS_build38 = as.integer(start(drive_gr38)))]
  }

  drive_variant_info <- merge(drive_variant_info, drive_liftover_map, by = "row_id", all.x = TRUE, sort = FALSE)[order(row_id)]

  drive_variant_info[, `:=`(
    CHR_build37 = CHR,
    POS_build37 = POS,
    ID_build37 = ID,
    ID_build38 = fifelse(n_lifted == 1L, paste(CHR_build38, POS_build38, REF, ALT, sep = ":"), NA_character_)
  )]

  drive_map_build38 <- drive_variant_info[
    n_lifted == 1L,
    .(
      row_id,
      CHR_build37,
      POS_build37,
      ID_build37,
      CHR_build38,
      POS_build38,
      ID_build38,
      REF,
      ALT,
      QUAL,
      FILTER,
      INFO,
      FORMAT
    )
  ]

  drive_map_unmapped <- drive_variant_info[
    n_lifted != 1L,
    .(
      row_id,
      CHR_build37,
      POS_build37,
      ID_build37,
      REF,
      ALT,
      n_lifted
    )
  ]

  fwrite(x = drive_map_build38, file = map_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(x = drive_map_unmapped, file = unmapped_file, sep = "\t", quote = FALSE, na = "NA")

  liftover_summary[[chr]] <- data.table(
    CHR = chr,
    input_variants = nrow(drive_variant_info),
    mapped_variants = nrow(drive_map_build38),
    unmapped_variants = nrow(drive_map_unmapped),
    mapping_rate = nrow(drive_map_build38) / nrow(drive_variant_info),
    status = "completed"
  )

  message("Chromosome ", chr, ": ", nrow(drive_map_build38), " of ", nrow(drive_variant_info), " positions mapped")

  rm(drive_variant_info, drive_liftover_map, drive_gr37, drive_gr38_list, drive_map_build38, drive_map_unmapped)
  if (exists("drive_gr38")) rm(drive_gr38)
  gc()
}

liftover_summary <- rbindlist(liftover_summary, use.names = TRUE, fill = TRUE)

fwrite(
  x = liftover_summary,
  file = file.path(drive_geno_dir, "oncoarray_build37_to_build38_summary.csv")
)

print(liftover_summary)
