#' Calculate Fst estimates for three sets of variants via PLINK
#' https://www.cog-genomics.org/plink2/basic_stats

#' Given a subset of sub-populations defined via --within, --fst writes FST
#' estimates for each autosomal diploid variant (computed using the method
#' introduced in Weir BS, Cockerham CC (1984) Estimating F-statistics for the
#' analysis of population structure) to plink.fst, and reports raw and weighted
#' global means to the log.

#' @param geno1F (char) path to file with SNP genotype data (PLINK format)
#' @param geno1_fam1F (char) path to file with SNP genotype data for the
#' control population (eg. CEU)
#' @param geno1_fam2F (char) path to file with SNP genotype data for the
#' case population (eg. ASW)
#' @param fstModifier (char) if you only have two sub-populations, you can
#' represent them with case/control status instead of clusters, and use the
#' case-control' modifier to request Fst estimates based on them (optional)
#' @param snpGroup (string) string indicating the 3 SNP groups that Fst
#' will be calculated for (e.g., 'top' GSEA pathway SNPs)
#' @return (char) PLINK tables containing Fst estimates
#' @export

calculateFst <- function(geno1F, geno1_fam1F, geno1_fam2F,
                         fstModifier, snpGroup=c("top","bottom","random"),
                         topOutDir, bottomOutDir, randomOutDir) {

  # Load SNP lists for each top pathway
  topPaths <- list.files(path=sprintf("%s/pathways_top", pathwaySnpDir),
                         pattern="*", full.names=T, recursive=F)

  for (i in 1:length(topPaths)) {

    # Subset original PLINK file per each SNP list
    if (!file.exists(topPlinkDir)) dir.create(topPlinkDir)
    str1 <- sprintf("%s --bfile %s --extract %s", PLINK, geno1F, topPaths[i])
    str2 <- sprintf("--make-bed --allow-no-sex --out %s/%s",
                    topPlinkDir, basename(topPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)

    # Calculate Fst for each set of top pathway SNPs
    if (!file.exists(topOutDir)) dir.create(topOutDir)
    str1 <- sprintf("%s --bfile %s/%s --fst %s", PLINK, topPlinkDir,
                    basename(topPaths[i]), fstModifier)
    str2 <- sprintf("--allow-no-sex --out %s/%s", topOutDir,
                    basename(topPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  # Load SNP lists for each bottom pathway
  bottomPaths <- list.files(path=sprintf("%s/pathways_bottom", pathwaySnpDir),
                            pattern="*", full.names=T, recursive=F)

  for (i in 1:length(bottomPaths)) {

    # Subset original PLINK file per each SNP list
    if (!file.exists(bottomPlinkDir)) dir.create(bottomPlinkDir)
    str1 <- sprintf("%s --bfile %s --extract %s", PLINK, geno1F,bottomPaths[i])
    str2 <- sprintf("--make-bed --allow-no-sex --out %s/%s",
                    bottomPlinkDir, basename(bottomPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)

    # Calculate Fst for each set of bottom pathway SNPs
    if (!file.exists(bottomOutDir)) dir.create(bottomOutDir)
    str1 <- sprintf("%s --bfile %s/%s --fst %s", PLINK, bottomPlinkDir,
                    basename(bottomPaths[i]), fstModifier)
    str2 <- sprintf("--allow-no-sex --out %s/%s", bottomOutDir,
                    basename(bottomPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  # Load SNP lists for each randomly selected pathway
  randomPaths <- list.files(path=sprintf("%s/pathways_random", pathwaySnpDir),
                            pattern="*", full.names=T, recursive=F)

  for (i in 1:length(randomPaths)) {

    # Subset original PLINK file per each SNP list
    if (!file.exists(randomPlinkDir)) dir.create(randomPlinkDir)
    str1 <- sprintf("%s --bfile %s --extract %s", PLINK, geno1F,randomPaths[i])
    str2 <- sprintf("--make-bed --allow-no-sex --out %s/%s",
                   randomPlinkDir, basename(randomPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)

    # Calculate Fst for each set of random pathway SNPs
    if (!file.exists(randomOutDir)) dir.create(randomOutDir)
    str1 <- sprintf("%s --bfile %s/%s --fst %s", PLINK, randomPlinkDir,
                   basename(randomPaths[i]), fstModifier)
    str2 <- sprintf("--allow-no-sex --out %s/%s", randomOutDir,
                    basename(randomPaths[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  ## Generate dummy data for negative controls using original PLINK data
  ## within-population Fst calculation (expecting low values since allele
  ## frequencies shouldn't vary within a sub-population)
  famString <- c(geno1_fam1F, geno1_fam2F)

  for (i in 1:length(famString)) {
    cat("Generating dummy data for negative controls...\n")
    fam <- read.table(file=sprintf("%s.fam", famString[i]), h=F, as.is=T)
    all_samples <- seq_len(nrow(fam))
    sample_1 <- sample(all_samples, length(all_samples)/2) #randomly select half of samples
    sample_2 <- all_samples[!all_samples %in% sample_1]
    split_1 <- fam[sample_1, ]; split_1[,6] <- 1 #code half '1'
    split_2 <- fam[sample_2, ]; split_2[,6] <- 2 #code other half '2'
    dummy_split <- rbind(split_1, split_2)
    table(dummy_split[,6])
    write.table(dummy_split, file=sprintf("%s/tmp.fam",
                dirname(famString[i])), col=F, row=F, quote=F, sep="\t")

    # Copy original PLINK .bim and .bed files to tmp files
    system(sprintf("cp %s.bim %s/tmp.bim", famString[i], dirname(famString[i])))
    system(sprintf("cp %s.bed %s/tmp.bed", famString[i], dirname(famString[i])))

    # Calculate within-population Fst
    str1 <- sprintf("%s --bfile %s/tmp --fst %s", PLINK, dirname(famString[i]),
                    fstModifier)
    str2 <- sprintf("--allow-no-sex --out %s/tmp", dirname(famString[i]))
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  # PLOT PATHWAY FST VALS
  # c.bind.fill() combines lists with uneven rows (work-around for cbind())
  ## For example, top pathway 1 has 83 rows (SNPs) and top pathway 2 has
  ## 111 rows (SNPs). cbind() isn't normally able to combine lists with
  ## an uneven number of rows, thus cbind.fill() will concatenate them
  ## without losing the original data structure
  #### nrow(top_fst[[1]]) -> 83    length(top_fst_list[[1]]) -> 83
  #### nrow(top_fst[[2]]) -> 111   length(top_fst_list[[2]]) -> 111

  cbind.fill <- function(...) {
    transpoted <- lapply(list(...),t)
    transpoted_dataframe <- lapply(transpoted, as.data.frame)
    return(data.frame(t(rbind.fill(transpoted_dataframe))))
  }

  ##top
  top_fst <- list.files(topOutDir, pattern="*.fst$", full.names=T)
  top_fst <- lapply(top_fst, function(x)read.table(x, h=T))
  top_fst_list <- c()
  for(i in 1:length(top_fst)) {
    top_fst_list[i] <- do.call(cbind.fill, list(top_fst[[i]]$FST))
  }

  ##bottom
  bottom_fst <- list.files(bottomOutDir, pattern="*.fst$", full.names=T)
  bottom_fst <- lapply(bottom_fst, function(x)read.table(x, h=T))
  bottom_fst_list <- c()
  for(i in 1:length(bottom_fst)) {
    bottom_fst_list[i] <- do.call(cbind.fill, list(bottom_fst[[i]]$FST))
  }

  ##random
  random_fst <- list.files(randomOutDir, pattern="*.fst$", full.names=T)
  random_fst <- lapply(random_fst, function(x)read.table(x, h=T))
  random_fst_list <- c()
  for(i in 1:length(random_fst)) {
    random_fst_list[i] <- do.call(cbind.fill, list(random_fst[[i]]$FST))
  }

  ##controls
  control_ceu <- read.table(sprintf("%s/tmp.fst", dirname(geno1_fam1F)),
                  h=T, as.is=T)
  control_asw <- read.table(sprintf("%s/tmp.fst", dirname(geno1_fam2F)),
                  h=T, as.is=T)
  ceu_fst_list <- list(control_ceu$FST)
  asw_fst_list <- list(control_asw$FST)

  all_fst_list <- c(ceu_fst_list, asw_fst_list, top_fst_list,
                    bottom_fst_list, random_fst_list)

  #top_path_names <- gsub("\\%.*", "", basename(topPaths)) #remove characters after '%'
  #bottom_path_names <- gsub("\\%.*", "", basename(bottomPaths))
  #random_path_names <- gsub("\\%.*", "", basename(randomPaths))

  cat("Generating plots...\n")
  pdf(file=sprintf("%s/fst_boxplots.pdf", outDir))
  op <- par(mar = c(4,4,2,1) + 0.1, cex.axis=0.8, cex.lab=0.9)
  #boxplot(ceu_fst_list, ylab="Fst", xlab="CEU", col="thistle1",
  #        main="Variation of SNPs within the European sub-population")
  #boxplot(asw_fst_list, ylab="Fst", xlab="ASW", col="thistle4",
  #        main="Variation of SNPs within the African sub-population")
  boxplot(top_fst_list, ylab="Fst", xlab="Pathways", col="mediumpurple1",
          main=paste0("Varation of SNPs in top pathways\n",
                      "between European and African populations"))
  abline(h=0.1, col="red")
  boxplot(bottom_fst_list, ylab="Fst", xlab="Pathways", col="mediumseagreen",
          main=paste0("Variation of SNPs in bottom pathways\n",
                      "between European and African populations"))
  abline(h=0.1, col="red")
  boxplot(random_fst_list, ylab="Fst", xlab="Pathways", col="lightsalmon",
          main=paste0("Variation of SNPs in randomly selected",
                      "pathways\n between European and African population"))
  abline(h=0.1, col="red")
  par(op) ##reset


  #######WEIGHTED MEAN##################3
  # Compare median of each SNP group
  ceu_fst_df <- as.data.frame(unlist(ceu_fst_list))
  asw_fst_df <- as.data.frame(unlist(asw_fst_list))
  top_fst_df <- as.data.frame(unlist(top_fst_list))
  bottom_fst_df <- as.data.frame(unlist(bottom_fst_list))
  random_fst_df <- as.data.frame(unlist(random_fst_list))
  all_fst_df <- c(ceu_fst_df, asw_fst_df, top_fst_df,
		              bottom_fst_df, random_fst_df)
  boxplot(all_fst_df, xlab="Fst", ylab="Group",
          main="Mean variation of SNPs among different population groups")
  dev.off()
  cat(" done.\n")

  #Clean up
  unlink(sprintf("%s/tmp.*", dirname(geno1_fam1F)))
  unlink(sprintf("%s/tmp.*", dirname(geno1_fam2F)))
}

#highuv <- which(top_data_list[[14]]$FST > 0.4)
#top_data_list[[14]] [highuv,]

#CHR       SNP       POS NMISS      FST
#55  11 rs4082224 113873697  4945 0.447754
