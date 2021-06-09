#' Script to calculate population stratification between different
#' ancestry groups (e.g., CEU vs ASW) via complete linkage clustering

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param famF (char) path to PLINK .fam file (case/control coded)
#' @param outF (char) path to save PCA .png image
#' @param dimensions (integer)
#' @export

popPCA <- function(genoF, famF, PLINK, dimensions=3L, outF) {

  # Calculate PCA via PLINK
#  str1 <- sprintf("%s --bfile %s --allow-no-sex", PLINK, genoF)
  str1 <- sprintf("%s --bed %s.bed --bim %s.bim --fam %s", PLINK, genoF,
                  genoF, famF)
  str2 <- sprintf("--cluster --mds-plot %i, --allow-no-sex --out %s", dimensions, outF)
  cmd <- paste(str1, str2)
  system(cmd)

  # Plot PCA
  ## Read in MDS (multidimensional scaling) table
  mds <- read.table(sprintf("%s.mds", outF), h=T, as.is=T)

  ## Read in subject information pertaining to the different ancestry groups
  ## (NOTE: PLINK .fam file)
  pops <- read.table(famF, h=F, as.is=T)
  pops.ctrls <- filter(pops, V6 == 1)
  pops.cases <- filter(pops, V6 == 2)

  ## Identify MDS components per ancestry groups for cluster grouping
  for (i in 1:nrow(mds)) {
    if (mds$IID[i] %in% pops.ctrls[,2]) mds[i, "population"] <- pop1
    else if (mds$IID[i] %in% pops.cases[,2]) mds[i, "population"] <- pop2
    else mds[i, "population"] <- "Other"
  }
  mds$population <- as.factor(mds$population)

  ## Plot results via ggplot
  ggplot(mds, aes(x=C1, y=C2, colour=population)) +
      geom_point(shape=4, aes(size=population)) +
      ggtitle(ggtitle) +
      labs(x="PC1", y="PC2") +
      ylim(-0.10, 0.15) +
      scale_color_manual(values=c("#fdb462", "#984ea3", "lightgrey")) +
      scale_size_manual(values=c(2,2,0.8)) +
      theme_Publication()

  ggsave(sprintf("%s.png", outF), width=7, height=7)

}
