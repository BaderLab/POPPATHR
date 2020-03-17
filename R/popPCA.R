#' Determines population stratification between different
#' ancestry groups (e.g., CEU vs YRI) via complete linkage clustering
#'
#' @param genotype_file (char) path to PLINK-formatted SNP genotype data.
#' @param fam_file (char) path to population-coded PLINK fam file.
#' @param pop_one (char) character code for the first population (e.g., CEU).
#' @param pop_two (char) character code for the second population (e.g., YRI).
#' @param dimensions (integer) dimensions for which to calculate PCA sim (default=3).
#' @param out_file (char) path to write PNG image.
#'
#' @return none
#' @export
#'

popPCA <- function(genotype_file, fam_file, pop_one, pop_two,
                   dimensions=3, out_file) {

  # Calculate PCA via PLINK
  str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s", genotype_file, genotype_file, fam_file)
  str2 <- sprintf("--cluster --mds-plot %i --allow-no-sex --out %s", dimensions, out_file)
  cmd <- paste(str1, str2)
  system(cmd)

  # Read in MDS (multidimensional scaling) table
  mds <- read.table(sprintf("%s.mds", out_file), h=TRUE, as.is=TRUE)

  # Read in subject information pertaining to the different ancestry groups
  pops <- read.table(fam_file, h=FALSE, as.is=TRUE)
  pops.ctrls <- filter(pops, V6 == 1)
  pops.cases <- filter(pops, V6 == 2)

  # Identify MDS components per ancestry groups for cluster grouping
  for (i in 1:nrow(mds)) {
    if (mds$IID[i] %in% pops.ctrls[,2]) mds[i, "population"] <- pop_one
    else if (mds$IID[i] %in% pops.cases[,2]) mds[i, "population"] <- pop_two
    else mds[i, "population"] <- "Other populations"
  }
  mds$population <- as.factor(mds$population)

  # Plot PCA
  ggplot(mds, aes(x=C1, y=C2, colour=population)) +
      geom_point(shape=4, size=2) +
      ggtitle("Population stratification via PCA") +
      labs(x="PC1", y="PC2") +
      scale_color_brewer(palette="Dark2") +
      theme_classic()
  ggsave(sprintf("%s.png", out_file), width=5, height=5)
}
