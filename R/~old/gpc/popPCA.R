#' Script to calculate population stratification between different
#' ancestry groups (e.g., CEU vs ASW) via complete linkage clustering

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param famF (char) path to PLINK .fam file (case/control coded)
#' @param outF (char) path to save PCA .png image
#' @param dimensions (integer)

popPCA <- function(genoF, famF, PLINK, dimensions=3L, outF) {

  # Calculate PCA via PLINK
  str1 <- sprintf("%s --bfile %s --allow-no-sex", PLINK, genoF)
  str2 <- sprintf("--cluster --mds-plot %i, --out %s", dimensions, outF)
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
    else mds[i, "population"] <- NA
  }

  ## Plot results via ggplot
  ggplot(mds, aes(x=C1, y=C2, colour=population)) +
      geom_point() +
      ggtitle(ggtitle) +
      scale_x_continuous(name="PC1") +
      scale_y_continuous(name="PC2") +
      scale_color_discrete(breaks=c(pop1, pop2)) + #keep only 'CEU' + 'ASW' in legend
      geom_hline(yintercept=0, colour="grey", linetype="dashed", size=1) +
      geom_vline(xintercept=0, colour="grey", linetype="dashed", size=1) +
      theme_set(theme_minimal()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank())

  ggsave(sprintf("%s.png", outF), width=9, height=9)

}
