#' Calculates SNP-level FST values given two populations
#'
#' @param genotype_file (char) path to PLINK-formatted SNP genotype data.
#' @param fam_file (char) path to PLINK population-coded fam file.
#' @param output_folder (char) path to output directory.
#' @param output_file (char) path to final GSEA-formatted input file.
#'
#' @return none
#' @export
#'

calcFST <- function(genotype_file, fam_file, output_folder, output_file) {

  # Set up PLINK to calculate FST between two populations
  cat(sprintf("* Using PLINK to calculate FST between %s and %s populations\n", pop_one, pop_two))
  str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s", genotype_file, genotype_file, fam_file)
  str2 <- sprintf("--fst case-control --allow-no-sex")
  str3 <- sprintf("--out %s", sprintf("%s/%s", output_folder, basename(genotype_file)))
  cmd  <- sprintf("%s %s %s", str1, str2, str3)
  system(cmd)

  # Read in generated file
  fst_file <- sprintf("%s/%s.fst", output_folder, basename(genotype_file))
  fst_cc <- fread(fst_file, h=TRUE, data.table=FALSE)
  fst_input <- subset(fst_cc, select=c(SNP, FST))
  fst_input[is.na(fst_input)] <- 0

  # NOTE: header of column 2 is labelled "CHI2" for purposes of reading into
  # calculate_gsea.pl. This script explicitly accepts two header labels,
  # 'P' or 'CHI2' typically from GWAS. The column still contains SNP FST values.
  cat(sprintf("* Writing out FST file for %s SNPs to %s\n", nrow(fst_cc), output_file))
  write.table(fst_input, file=output_file, col.names=c("Marker", "CHI2"), row=FALSE, sep="\t", quote=FALSE)
}
