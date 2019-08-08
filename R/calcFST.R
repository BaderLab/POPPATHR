#' Calculates SNP-level FST values for use in runGSEA.R
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format).
#' @param realFam (char) path to PLINK case/control coded fam file.
#' @param outDir (char) directory to store output files.
#' @param outF (char) path to write SNP-FST file.
#'
#' @return none
#' @export
#'
calcFST <- function(genoF, realFam, outDir, outF) {
  # Set up PLINK to calculate MAF between case/control populations
  str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s", genoF, genoF, realFam)
  str2 <- sprintf("--fst case-control --allow-no-sex")
  str3 <- sprintf("--out %s", sprintf("%s/%s", outDir, basename(genoF)))
  cmd  <- sprintf("%s %s %s", str1, str2, str3)
  system(cmd)

  # Read in generated plink.fst file
  fst_cc  <- fread(sprintf("%s/%s.fst", outDir, basename(genoF)),
		h=TRUE, data.table=FALSE)
  fst_input <- subset(fst_cc, select=c(SNP, FST))
  fst_input[is.na(fst_input)] <- 0

  # NOTE: header of column 2 is labelled "CHI2" for purposes of reading into
  # calculate_gsea.pl. This script explicitly accepts two header labels,
  # 'P' or 'CHI2' typically from GWAS. The column still contains SNP FST values.
  cat(sprintf("*Writing out FST file to %s.\n", outF))
  write.table(fst_input, file=outF,
		col.names=c("Marker", "CHI2"), row=FALSE, sep="\t", quote=FALSE)
}
