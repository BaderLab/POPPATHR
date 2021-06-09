#' Calculate R2 estimates for three sets of variants via PLINK
#' https://www.cog-genomics.org/plink2/basic_stats

#' @param geno1F (char) path to file with SNP genotype data (PLINK format)
#' @param pathwaySnpDir (char) path to directory with pathway SNP lists
#' @param LDwindow (integer) value for PLINK --ld-window-r2
#' @return (char) PLINK table containing Fst estimates
#' @export

calculateR2 <- function(geno1F, pathwaySnpDir, LDwindow=0L,
                        topOutDir, bottomOutDir, randomOutDir) {

  topPaths <- list.files(path=sprintf("%s/pathways_top", pathwaySnpDir),
                         pattern="*", full.names=F, recursive=F)

  for (i in 1:length(topPaths)) {
    if (!file.exists(topOutDir)) dir.create(topOutDir)
    str1 <- sprintf("%s --bfile %s/%s --r2 inter-chr --ld-window-r2 %i",
                    PLINK, topPlinkDir, topPaths[i], LDwindow)
    str2 <- sprintf("--allow-no-sex --out %s/%s", topOutDir, topPaths[i])
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  bottomPaths <- list.files(path=sprintf("%s/pathways_bottom", pathwaySnpDir),
                         pattern="*", full.names=F, recursive=F)

  for (i in 1:length(bottomPaths)) {
    if (!file.exists(bottomOutDir)) dir.create(bottomOutDir)
    str1 <- sprintf("%s --bfile %s/%s --r2 inter-chr --ld-window-r2 %i",
                    PLINK, bottomPlinkDir, bottomPaths[i], LDwindow)
    str2 <- sprintf("--allow-no-sex --out %s/%s", bottomOutDir, bottomPaths[i])
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

  randomPaths <- list.files(path=sprintf("%s/pathways_random", pathwaySnpDir),
                            pattern="*", full.names=F, recursive=F)

  for (i in 1:length(randomPaths)) {
    if (!file.exists(randomOutDir)) dir.create(randomOutDir)
    str1 <- sprintf("%s --bfile %s/%s --r2 inter-chr --ld-window-r2 %i",
                    PLINK, randomPlinkDir, randomPaths[i], LDwindow)
    str2 <- sprintf("--allow-no-sex --out %s/%s", randomOutDir, randomPaths[i])
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)
  }

    # PLOT PATHWAY R2 VALS
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
      return (data.frame(t(rbind.fill(transpoted_dataframe))))
    }

    ##top
    top_r2 <- list.files(topOutDir, pattern="*.ld", full.names=T)
    top_r2 <- lapply(top_r2, function(x)read.table(x, h=T))
    top_r2_list <- c()
    for(i in 1:length(top_r2)) {
      top_r2_list[i] <- do.call(cbind.fill, list(top_r2[[i]]$R2))
    }

    ##bottom
    bottom_r2 <- list.files(bottomOutDir, pattern="*.ld", full.names=T)
    bottom_r2 <- lapply(bottom_r2, function(x)read.table(x, h=T))
    bottom_r2_list <- c()
    for(i in 1:length(bottom_r2)) {
      bottom_r2_list[i] <- do.call(cbind.fill, list(bottom_r2[[i]]$R2))
    }

    ##random
    random_r2 <- list.files(randomOutDir, pattern="*.ld", full.names=T)
    random_r2 <- lapply(random_r2, function(x)read.table(x, h=T))
    random_r2_list <- c()
    for(i in 1:length(random_r2)) {
      random_r2_list[i] <- do.call(cbind.fill, list(random_r2[[i]]$R2))
    }

    #top_path_names <- gsub("\\%.*", "", topPaths) #remove characters after '%'
    #bottom_path_names <- gsub("\\%.*", "", bottomPaths)
    #random_path_names <- gsub("\\%.*", "", randomPaths)

    # Save separate plots within pdf doc
    cat("Generating plots...\n")
    pdf(file=sprintf("%s/r2_boxplots.pdf", outDir))
    op <- par(mar = c(4,4,2,1) + 0.1, cex.axis=0.8, cex.lab=0.9)
    boxplot(top_r2_list, ylab="R2", xlab="Pathways", col="mediumpurple1",
            main="R2 distribution for SNPs in top pathways")
    abline(h=0.2, col="red")
    boxplot(bottom_r2_list, ylab="R2", xlab="Pathways", col="mediumseagreen",
            main="R2 distribution for SNPs in bottom pathways")
    abline(h=0.2, col="red")
    boxplot(random_r2_list, ylab="R2", xlab="Pathways", col="lightsalmon",
            main="R2 distribution for SNPs in randomly selected pathways")
    abline(h=0.2, col="red")
    par(op) ##reset

    # Filter for SNP pairs with R2 > 0.2
    cat("Filtering SNP pairs at R2 = 0.2\n")
    top_r2_df <- as.data.frame(unlist(top_r2_list))
    top_r2_df_filter <- subset(top_r2_df, unlist(top_r2_list) > 0.2)
    cat(sprintf("Top: %i from %i pairs\n", nrow(top_r2_df_filter),
                nrow(top_r2_df)))

    bottom_r2_df <- as.data.frame(unlist(bottom_r2_list))
    bottom_r2_df_filter <- subset(bottom_r2_df, unlist(bottom_r2_list) > 0.2)
    cat(sprintf("Bottom: %i from %i pairs\n", nrow(bottom_r2_df_filter),
                nrow(bottom_r2_df)))

    random_r2_df <- as.data.frame(unlist(random_r2_list))
    random_r2_df_filter <- subset(random_r2_df, unlist(random_r2_list) > 0.2)
    cat(sprintf("Random: %i from %i pairs\n", nrow(random_r2_df_filter),
                nrow(random_r2_df)))

    # Plot median R2 distribution per group
    all_r2_df_filter <- c(top_r2_df_filter, bottom_r2_df_filter,
                          random_r2_df_filter)
    boxplot(all_r2_df_filter, ylab="R2", xlab="Group",
            col=c("mediumpurple1", "mediumseagreen", "lightsalmon"),
            main="Median R2 (cutoff=0.2) for SNPs across sample groups")
    violins(all_r2_df_filter, connect = c("median", "mean"), ylab="R2",
            xlab="Group", names=1:3,
            col=c("mediumpurple1", "mediumseagreen", "lightsalmon"),
            main="Median R2 (cutoff=0.2) for SNPs across sample groups")
    dev.off()
    cat(" done.\n")
}
