#' Calculates difference in minor allele frequency between two populations
#' for use in runGSEA.R

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param fam1F (char) path to PLINK .fam file for population 1
#' @param fam2F (char) path to PLINK .fam file for population 2
#' @return (char) data frame with two columns: 1) SNP ID (marker),
#' 2) SNP association statistic
#' @export

calculateMAFdiff <- function(genoF, fam1F, fam2F) {
out <- tryCatch({

  # Setup plink to calculate MAF for population 1
  cat("Calculating minor allele frequency for population 1 via PLINK...")
  str1 <- sprintf("%s --bfile %s ", PLINK, genoF)
  str2 <- sprintf("--freq --allow-no-sex --within %s", fam1F)
  str3 <- sprintf("--out %s/%s", plinkDir, pop1)
  cmd  <- sprintf("%s %s %s", str1, str2, str3)
  system(cmd)
  cat(sprintf(" log file written to %s/%s.log\n", plinkDir, pop1))

  # Setup plink to calculate MAF for population 2
  cat("Calculating minor allele frequency for population 2 via PLINK...")
  str1 <- sprintf("%s --bfile %s", PLINK, genoF)
  str2 <- sprintf("--freq --allow-no-sex --within %s", fam2F)
  str3 <- sprintf("--out %s/%s", plinkDir, pop2)
  cmd  <- sprintf("%s %s %s", str1, str2, str3)
  system(cmd)
  cat(sprintf(" log file written to %s/%s.log\n", plinkDir, pop2))

  # Reformat frequency files for HM3 dataset
  ### Individuals in the HM3 dataset are coded according to their Family ID,
  ### Paternal ID, and Maternal ID. Running PLINK --freq along with the
  ### --within option allows for any categorical grouping of the sample.
  ### In the previous step, we are clustering the samples via their ethnicity.
  ### However, since the HM3 variants are also coded via 3 distinct IDs,
  ### PLINK additionally categorizes the variants based on those IDs

  ### head -n5 [pop].frq.strat
  ### CHR           SNP     CLST   A1   A2      MAF    MAC  NCHROBS
  ###   1     rs4124251        0    C    T        0      0        0
  ###   1     rs4124251  NA19700    C    T        0      0        0
  ###   1     rs4124251  NA19703    C    T        0      0        0
  ###   1     rs4124251  NA19818    C    T        0      0        0

  ### After removing duplicate SNP MAF calculations due to ID clustering
  ### head -n5 [pop].frq.strat2
  ### CHR          SNP     CLST   A1   A2      MAF    MAC  NCHROBS
  ### 1      rs4124251        0    C    T        0      0        0
  ### 1      rs6650104        0    G    A        0      0        0
  ### 1     rs10458597        0    T    C        0      0        0
  ### 1      rs9629043        0    T    C  0.02041      2       98

  ### In this manner, we effectively re-format the .frq.strat file
  ### to give us a single MAF calculation per SNP

  cat("Formatting PLINK frequency files...")
  pop1_frq  <- system(sprintf(paste("awk '(NR == 1) || ($3 == 0)'",
                                    "%s/%s.frq.strat > %s/%s.frq.strat2"),
                                    plinkDir, pop1, plinkDir, pop1))
  pop2_frq  <- system(sprintf(paste("awk '(NR == 1) || ($3 == 0)'",
                                    "%s/%s.frq.strat > %s/%s.frq.strat2"),
                                    plinkDir, pop2, plinkDir, pop2))
  cat(sprintf(" formatted files written to %s/*.frq.strat2\n", plinkDir))

  # Read in generated plink frequency files and modify column names
  pop1_frq  <- read.table(sprintf("%s/%s.frq.strat2", plinkDir, pop1),
                                  h=T, as.is=T)
  colnames(pop1_frq)[3:8] <- paste(colnames(pop1_frq[,c(3:8)]),
                                            "pop1", sep="_")

  pop2_frq  <- read.table(sprintf("%s/%s.frq.strat2", plinkDir, pop2),
                                  h=T, as.is=T)
  colnames(pop2_frq)[3:8] <- paste(colnames(pop2_frq[,c(3:8)]),
                                            "pop2", sep="_")

  # Check if the minor allele in both populations are the same!
  cat("Checking if minor alleles match between both populations...")
  pop1_pop2_frq  <- merge(pop1_frq, pop2_frq, by=c("SNP", "CHR"))
  check_same_MA <- all(pop1_pop2_frq$A1_pop1 == pop1_pop2_frq$A1_pop2)
  cat(sprintf(" %s.\n", check_same_MA))

# Stop program if minor alleles don't match
  if (check_same_MA != TRUE) {

      stop("Warning: minor alleles in both populations don't match!\n")

  } else {

  # Calculate the difference bewteen compared population MAFs
  both_MAF   <- subset(pop1_pop2_frq, select=c(SNP, CHR, MAF_pop1, MAF_pop2))
  ## ABSOLUTE MAF ##
  both_MAF$MAFdiff <- both_MAF$MAF_pop1 - both_MAF$MAF_pop2
  marker_MAF_input  <- subset(both_MAF, select=c(SNP, MAFdiff))

  }

  # NOTE: header of column 2 is labelled "CHI2" for purposes of reading into
  # calculate_gsea.pl. This script explicitly accepts two header labels,
  # 'P' or 'CHI2'. The column still contains MAF difference values
  cat(sprintf("Writing out GSEA input file to %s/markerMAF.txt\n", plinkDir))
  write.table(marker_MAF_input, file=sprintf("%s/markerMAF.txt", plinkDir),
              col.names=c("Marker", "CHI2"), row=F, sep="\t", quote=F)

  })
  return(out)
  cat("...closing log.\n")
}
