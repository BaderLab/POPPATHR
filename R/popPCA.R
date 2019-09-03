#' Calculate population stratification between different
#' ancestry groups (e.g., CEU vs YRI) via complete linkage clustering
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format).
#' @param famF (char) path to PLINK fam file (case/control population coded).
#' @param pop1 (char) character code for the first population (controls).
#' @param pop2 (char) character code for the second population (cases).
#' @param dimensions (integer, default=3).
#' @param outF (char) path to write PNG image.
#'
#' @return none
#' @import dplyr, ggplot2
#' @export
#'

popPCA <- function(genoF, famF, pop1, pop2, dimensions=3L, outF) {
  # Calculate PCA via PLINK
  str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s", genoF, genoF, famF)
  str2 <- sprintf("--cluster --mds-plot %i --allow-no-sex --out %s", dimensions, outF)
  cmd <- paste(str1, str2)
  system(cmd)

  # Plot PCA
  ## Read in MDS (multidimensional scaling) table
  mds <- read.table(sprintf("%s.mds", outF), h=TRUE, as.is=TRUE)

  ## Read in subject information pertaining to the different ancestry groups
  ## (NOTE: PLINK .fam file)
  pops <- read.table(famF, h=FALSE, as.is=TRUE)
  pops.ctrls <- filter(pops, V6 == 1)
  pops.cases <- filter(pops, V6 == 2)

  ## Identify MDS components per ancestry groups for cluster grouping
  for (i in 1:nrow(mds)) {
    if (mds$IID[i] %in% pops.ctrls[,2]) mds[i, "population"] <- pop1
    else if (mds$IID[i] %in% pops.cases[,2]) mds[i, "population"] <- pop2
    else mds[i, "population"] <- "Other populations"
  }
  mds$population <- as.factor(mds$population)

  ## Plot results via ggplot
  ggplot(mds, aes(x=C1, y=C2, colour=population)) +
      geom_point(shape=4, size=2) +
      ggtitle("Population stratification via PCA") +
      labs(x="PC1", y="PC2") +
      #ylim(-0.10, 0.15) +
      scale_color_brewer(palette="Dark2") +
      theme_classic()
  ggsave(sprintf("%s.png", outF), width=5, height=5)
}
