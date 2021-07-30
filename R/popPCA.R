#' Determines population stratification between different
#' ancestry groups (e.g., CEU vs YRI) via complete linkage clustering
#'
#' @param genotypeFile (char) path to PLINK-formatted SNP genotype data.
#' @param famFile (char) path to population-coded PLINK fam file.
#' @param popOne (char) character code for the first population (e.g., CEU).
#' @param popTwo (char) character code for the second population (e.g., YRI).
#' @param dimensions (integer) dimensions for which to calculate PCA sim (default=3).
#' @param outputFile (char) path to write PNG image.
#' @export
popPCA <- function(genotypeFile, famFile, popOne, popTwo,
    dimensions=3, outputFile) {

    # Use PLINK for principle component analysis (PCA)
    cat(sprintf("* Using PLINK to run PCA on %s and %s genotypes\n",
        popOne, popTwo))
    str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s",
        genotypeFile, genotypeFile, famFile)
    str2 <- sprintf("--cluster --mds-plot %i --allow-no-sex --out %s",
        dimensions, outputFile)
    cmd <- paste(str1, str2)
    system(cmd)

    # Read in MDS (multidimensional scaling) table
    mds <- read.table(sprintf("%s.mds", outputFile), header=TRUE, as.is=TRUE)

    # Read in subject information pertaining to the different ancestry groups
    populations <- read.table(famFile, header=FALSE, as.is=TRUE)
    populationsOne <- filter(populations, V6 == 1)
    populationsTwo <- filter(populations, V6 == 2)

    # Identify MDS components per ancestry groups for cluster grouping
    for (i in 1:nrow(mds)) {
        if (mds$IID[i] %in% populationsOne[,2]) {
            mds[i, "population"] <- popOne
        } else if (mds$IID[i] %in% populationsTwo[,2]) {
            mds[i, "population"] <- popTwo
        } else {
            mds[i, "population"] <- "Other populations"
        }
    }

    # Set factor level for plotting
    mds$population <- as.factor(mds$population)
    cols <- brewer.pal(3, "Dark2")
    cols[3] <- "grey80"

    # Plot PCA
    ggplot(mds, aes(x=C1, y=C2, colour=population)) +
        geom_point(shape=4, size=2) +
        ggtitle("Population stratification via PCA") +
        labs(x="PC1", y="PC2") +
        scale_color_manual(values=cols) +
        theme_classic(base_size=8)

    # Save plot
    cat(sprintf("* Drawing out PCA plot to %s.pdf\n", outputFile))
    ggsave(sprintf("%s.pdf", outputFile), width=4, height=3)
}
