#' Calculates SNP-level FST values for use in runGSEA.R
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param realFAM (char) path to PLINK case/control coded fam file
#' @param phenoPerm (logical) runs GSEA using phenotype permutations
#'		(i.e., randomly re-shuffing case/control labels n times).
#' 		if FALSE, will run GSEA using genotype permutations (see setupGSEArun.R)
#' @param permFAM (char) path to files with fam files from re-shuffled
#'    case/control labels. If phenoPerm is FALSE, leave this blank.
#' @param PLINK (char) path to PLINK executable
#' @param outDir (char) directory to store output files (PLINK frequency
#'    calculations plus GSEA input file)
#' @return (char) vector of data frame with two columns: 1) SNP ID (marker),
#'    2) SNP association statistic
#' @export

calcFST <- function(genoF, realFAM, PLINK,
                    phenoPerm, permFAM, outDir) {

  doCalc <- function(famF, plinkOut, frqOut) {

    # Setup PLINK to calculate MAF between case/control populations
    cat("* Calculating population-based FST\n")
    str1 <- sprintf("%s --bed %s.bed --bim %s.bim --fam %s", PLINK, genoF,
                    genoF, famF)
    str2 <- sprintf("--fst case-control --allow-no-sex")
    str3 <- sprintf("--out %s", plinkOut)
    cmd  <- sprintf("%s %s %s", str1, str2, str3)
    system(cmd)
    cat(" done.\n")

    # Read in generated plink.fst file and
    fst_cc  <- fread(sprintf("%s.fst", plinkOut), h=T, data.table=F)
    fst_input <- subset(fst_cc, select=c(SNP, FST))
    fst_input[is.na(fst_input)] <- 0

    # NOTE: header of column 2 is labelled "CHI2" for purposes of reading into
    # calculate_gsea.pl. This script explicitly accepts two header labels,
    # 'P' or 'CHI2'. The column still contains ΔMAF values
    cat(sprintf("Writing out FST file to %s.\n", frqOut))
    write.table(fst_input, file=frqOut,
                col.names=c("Marker", "CHI2"), row=F, sep="\t", quote=F)

  }

  # Real calculation
  doCalc(famF=realFAM,
         plinkOut=sprintf("%s/%s", outDir, basename(genoF)),
         frqOut=sprintf("%s/markerFST.txt", outDir))

  if (phenoPerm == TRUE) {
   # Perm calculations (shuffled case/control labels)
   # NOTE: using file_path_sans_ext() to ensure [k].fam corresponds with
   # [k].frq.cc and [k].txt (PLINK numerical order different from R)
   # i.e., with PLINK order is 100, 10, 11, ..., 19, 1, 20, 21, 22... etc
   # while with R length() order is 1, 2, 3, 4, 5, 6... etc
   perm_files <- list.files(path=permFAM, pattern="*.fam$", full.names=T)
   for (k in 1:length(perm_files)) {
     doCalc(famF=perm_files[k],
            plinkOut=file_path_sans_ext(perm_files[k]),
            frqOut=sprintf("%s.txt", file_path_sans_ext(perm_files[k])))
   }
 } else {
   return(NULL)
 }
}
