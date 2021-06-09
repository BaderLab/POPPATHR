#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param realFAM (char) path to PLINK case/control coded fam file
#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @param outDir (char) path to output directory
#' @return
#' @export

LDstats <- function(genoF, realFAM, hcInDir,
                    makePlots=TRUE, outDir) {

  # Calculate linkage disequilbrium statistics
  # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
  # to be considered, so that depth=1 would specify calculation of LD b/w
  # immediately adjacent SNPs. For our purposes we want to determine LD
  # b/w all SNPs in each pathway despite their distance from each other,
  # so we specify depth in ld() as ((number of SNPs)-1)
  calcLD <- function(snpDir, rep.num) {

    ss <- c(); snp.map <- c()
    ld.calc <- list(); pairwise.df <- list()
    diff.r2 <- list(); diff.dp <- list()

    bed <- sprintf("%s.bed", genoF)
    bim <- sprintf("%s.bim", genoF)
    fam <- realFAM

    cat(sprintf("Reading input SNP set %s... ", basename(genoF)))
    start.time <- Sys.time()
    all.ss <- read.plink(bed, bim, fam)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)

    ############################
    # experimenting with pop-based ld calc
  #  pop <- which(all.ss$fam$affected == 1)
  #  all.ss$genotypes <- all.ss$genotypes[pop, ]
    ############################

    # Get n SNPs per high-confidence pathways - gives the number to
    # sample for each N random set (N = # of hc pathways)
    snps <- list.files(path=snpDir, pattern="*.snps$", full.names=T)
    hc.paths <- lapply(snps, readLines)
    hc.path.size <- sapply(hc.paths, length)

#    rep.num <- rep.num
#    samples <- replicate(rep.num, hc.path.size)

    set.seed(1)
    for (i in 1:ncol(samples)) {
      n <- hc.path.size[i]

      cat(sprintf("random set %i with %i snps\n", i, n))
      ss[[i]] <- all.ss$genotypes[, sample(ncol(all.ss$genotypes),
                                  n, replace=F)]
      ld.calc[[i]] <- ld(ss[[i]],
                         stats=c("D.prime", "R.squared"),
                         depth=ncol(ss[[i]])-1)
      snp.map <- all.ss$map

      # Turn each LD matrix into a data frame
      r2 <- as.matrix(ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
      r2 <- subset(melt(r2), value!=0) #...for all non-zero values
      colnames(r2)[3] <- "R.squared"

      # Create dataframe containing pairwise distance calculations for each
      # SNP-SNP pair
      dp <- as.matrix(ld.calc[[i]]$D.prime)
      dp <- subset(melt(dp), value!=0)
      colnames(dp)[3] <- "D.prime"

      # Combine R2 and Dprime stats for each SNP-SNP pair
      all.stats <- merge(r2, dp, by=c("Var1", "Var2"))

      # Generate pariwise distance table for each SNP-SNP pair
      colnames(all.stats)[1] <- "snp.name"
      snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

      pairwise <- merge(snp.map, all.stats, by="snp.name")
      colnames(pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
      pairwise <- merge(snp.map, pairwise, by="snp.name")
      colnames(pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                              "chr_2", "pos_2", "R.squared", "D.prime")
      pairwise$dist <- abs(pairwise$pos_1 - pairwise$pos_2)
      pairwise$pathway <- sprintf("pathway_%i", i)

      pairwise.df[[i]] <- pairwise %>% mutate(R.squared=round(R.squared, 3))
      pairwise.df[[i]] <- pairwise %>% mutate(D.prime=round(D.prime, 3))

      diff.r2[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2) %>%
                             dplyr::select(R.squared) %>% unlist
      diff.dp[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2) %>%
                             dplyr::select(D.prime) %>% unlist
     }

    all.pairs <<- do.call("rbind", pairwise.df)
    diff.num <<- sapply(diff.r2, length)   #all SNP-SNP pairs per path
    diff.r2.mean <<- sapply(diff.r2, mean) #mean r2 value per pathway
    diff.dp.mean <<- sapply(diff.dp, mean) #mean dprime per pathway

    cat(sprintf("\nCalculated LD for %i total SNP pairs.\n", nrow(all.pairs)))
    cat(sprintf("\t --> %i total interchromosomal pairs.\n", sum(diff.num)))
  }

  cat("\nCalculating long-range correlation b/w randomly selected SNP sets.\n")
  calcLD(snpDir=hcInDir, rep.num=100L)
  rn.all.pairs <- all.pairs
  rn.diff.num <- diff.num
  rn.diff.r2.mean <- diff.r2.mean
  rn.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  #============================================================================#
  ## PLOT STATS
  if (makePlots == TRUE) {
    message("\nGenerating LD statistics plots.\n")

    # Set variables and other functions
    title <- paste("Degree of co-selection per interchromosomal SNP-SNP",
                   "\ninteraction within the high-confidence pathways",
                   "\nvs. randomly selected SNP groups")
    r2.axis.title <- bquote("Pairwise LD value (mean " *r^2*" of each SNP set)")
    dp.axis.title <- bquote("Pairwise LD value (mean D' of each SNP set)")

    ## for use with stat_summary(fun.data=box.style); allows white median line
    ## to appear after colouring and filling boxplots
    box.style <- function(x){
        return(c(y=median(x), ymin=median(x), ymax=median(x)))
    }

    ## for use with stat.summary(fun.data=give.n); displays sample size (N)
    ## courtesy of Bangyou at Stack Overflow
    give.n <- function(x){
        return(c(y=median(x)*1.50, label=length(x)))
        # experiment with the multiplier to find the perfect position
    }

    ## calculate empirical p by quantifying all permuted mean r2 values greater
    ## than the real mean r2 values, divided by the total number of replicates
    calc.p <- function(rn, hc){
      return( (length(which(rn > mean(hc)))+1) / (length(rn)+1) )
    }

    # Plotting mean R2 value for random SNP-SNP pairs on diff chromosome
    # against real mean R2 value for high-conf SNPs
    # interchromosomal SNP pairs used as proxy for co-selection/genetic ixns
    rn <- data.frame(R.squared=rn.diff.r2.mean,
                     D.prime=rn.diff.dp.mean,
                     pathway.group="rand")

    dist.dat <- rbind(rn)
    mean.dat <- ddply(dist.dat, "pathway.group", summarise,
                      R.squared.mean=mean(R.squared),
                      D.prime.mean=mean(D.prime))

    ggplot(dist.dat, aes(x=R.squared, colour=pathway.group,
                         fill=pathway.group)) +
        geom_density(alpha=0.3) +
        #geom_histogram(bins=20, alpha=0.5, position="identity") +
        geom_vline(data=mean.dat, aes(xintercept=R.squared.mean,
                                      colour=pathway.group),
                   linetype="dashed", size=1) +
        scale_x_continuous(r2.axis.title) +
        scale_y_continuous("Density") +
        ggtitle(title)
    ggsave(sprintf("%s/dist_hc-lc-rn_r2.png", outDir),
        width=8, height=7)

    ##with facets
    ggplot(dist.dat, aes(x=R.squared)) +
        facet_grid(pathway.group ~ .) +
        geom_density(alpha=0.3) +
        #geom_histogram(bins=20, alpha=0.5, position="identity") +
        geom_vline(data=mean.dat, aes(xintercept=R.squared.mean,
                                      colour=pathway.group),
                   linetype="dashed", size=1) +
        scale_x_continuous(r2.axis.title) +
        scale_y_continuous("Density") +
        ggtitle(title)
    ggsave(sprintf("%s/dist_hc-lc-rn_r2_facet.png", outDir),
        width=8, height=7)

    perm.p = calc.p(lc.diff.r2.mean, hc.diff.r2.mean)

    cat("\n-----------------------------------------------------------------")
    cat(sprintf("\nP-value for permuted sample vs. real test statistic = %g\n",
                perm.p))
       if (perm.p < 0.05) {
         cat("\nCelebrate! It's significant.")
       } else { cat("\nThat doesn't look so good.") }
    cat("\n-----------------------------------------------------------------\n")

    # Boxplot of total LD stat distribution b/w null vs. real
    ggplot(dist.dat, aes(x=pathway.group, y=R.squared)) +
          geom_boxplot(outlier.colour=NULL,
                       aes(colour=pathway.group, fill=pathway.group)) +
          stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
                       fun.data=box.style) +
          scale_y_continuous(r2.axis.title) +
          ggtitle(title)
          scale_x_discrete(labels=paste("N=", table(dist.dat$pathway.group),
                                        sep=""))
      ggsave(sprintf("%s/boxplot_hc-lc-rn_r2.png", outDir),
          width=8, height=7.5)
  }
 else { return() }
}
