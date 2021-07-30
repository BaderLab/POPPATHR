#' Calculates SNP-level FST values given two populations
#'
#' @param genotypeFile (char) path to PLINK-formatted SNP genotype data.
#' @param famFile (char) path to PLINK population-coded fam file.
#' @param popOne (char) character code for the first population (e.g., CEU).
#' @param popTwo (char) character code for the second population (e.g., YRI).
#' @param outputFolder (char) path to output directory.
#' @param outputFile (char) path to final GSEA-formatted input file.
#' @export
calcFST <- function(genotypeFile, famFile, popOne, popTwo,
    outputFolder, outputFile) {
    # Set up PLINK to calculate FST between two populations
    cat(sprintf("* Using PLINK to calculate FST between %s and %s populations\n",
        popOne, popTwo))
    str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s",
        genotypeFile, genotypeFile, famFile)
    str2 <- sprintf("--fst case-control --allow-no-sex")
    str3 <- sprintf("--out %s", sprintf("%s/%s", outputFolder,
    basename(genotypeFile)))
    cmd  <- sprintf("%s %s %s", str1, str2, str3)
    system(cmd)

    # Read in generated file
    fstFile <- sprintf("%s/%s.fst", outputFolder, basename(genotypeFile))
    fstCC <- fread(fstFile, header=TRUE, data.table=FALSE)
    fstInput <- subset(fstCC, select=c(SNP, FST))
    fstInput[is.na(fstInput)] <- 0

    # Write out file
    cat(sprintf("* Writing out FST file for %s SNPs to %s\n",
        nrow(fstCC), outputFile))
    write.table(fstInput, file=outputFile, col.names=c("snp", "fst"),
        row=FALSE, sep="\t", quote=FALSE)
}
