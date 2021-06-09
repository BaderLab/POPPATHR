#' Calculates difference in minor allele frequency between two populations
#' for use in runGSEA.R

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param fam1F (char) path to PLINK .fam file for population 1
#' @param fam2F (char) path to PLINK .fam file for population 2
#' @return (char) data frame with two columns: 1) SNP ID (marker),
#' 2) SNP association statistic
#' @export

calculateMAFdiff <- function(genoF) {

  # Setup plink to calculate MAF for population 1
  cat("Calculating case/control phenotype-stratified minor allele frequency...")
  str1 <- sprintf("%s --bfile %s ", PLINK, genoF)
  str2 <- sprintf("--freq case-control --allow-no-sex")
  str3 <- sprintf("--out %s/%s", plinkDir, basename(genoF))
  cmd  <- sprintf("%s %s %s", str1, str2, str3)
  system(cmd)
  cat(sprintf(" log file written to %s/%s.log\n", plinkDir, basename(genoF)))

  # Read in generated frequency file and calculate difference in minor
  # allele frequency for each SNP between both populations
  frq_cc  <- read.table(sprintf("%s/%s.frq.cc", plinkDir, basename(genoF)),
                                 h=T, as.is=T)

  ## ABSOLUTE MAF ##
  frq_cc$MAFdiff <- frq_cc$MAF_U - frq_cc$MAF_A
  marker_MAF_input  <- subset(frq_cc, select=c(SNP, MAFdiff))

  # NOTE: header of column 2 is labelled "CHI2" for purposes of reading into
  # calculate_gsea.pl. This script explicitly accepts two header labels,
  # 'P' or 'CHI2'. The column still contains MAF difference values
  cat(sprintf("Writing out GSEA input file to %s/markerMAF.txt\n", plinkDir))
  write.table(marker_MAF_input,
              file=sprintf("%s/markerMAF.txt", plinkDir),
              col.names=c("Marker", "CHI2"), row=F, sep="\t", quote=F)

}
