#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param lcInDir (char) path to files with low-confidence pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @param outDir (char) path to output directory
#' @return
#' @export

LDstats <- function(genoF, hcInDir, lcInDir,
                    makePlots=TRUE, outDir) {

  # Calculate linkage disequilbrium statistics
  # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
  # to be considered, so that depth=1 would specify calculation of LD b/w
  # immediately adjacent SNPs. For our purposes we want to determine LD
  # b/w all SNPs in each pathway despite their distance from each other,
  # so we specify depth in ld() as ((number of SNPs)-1)
  calcLD <- function(pathSet, snpDir) {

  ss <- c(); snp.map <- c()
  ld.calc <- list(); pairwise.df <- list()
  diff.r2 <- list(); diff.dp <- list()

    if (pathSet != "rand") {
      # Read PLINK files for each pathway
      bed <- list.files(path=snpDir, pattern="*.bed", full.names=T)
      bim <- list.files(path=snpDir, pattern="*.bim", full.names=T)
      fam <- list.files(path=snpDir, pattern="*.fam", full.names=T)

      # Convert PLINK files to snpStats input format
      # Output object is a list with 3 elements ($genotypes, $fam, $map)
      # NOTE: order is important!
      for (i in 1:length(bed)) {
        ss[[i]] <- read.plink(bed[i], bim[i], fam[i])

        ############################
        # experimenting with pop-based ld calc
        pop <- which(ss[[i]]$fam$affected == 2)
        ss[[i]]$genotypes <- ss[[i]]$genotypes[pop, ]
        ############################

        cat(sprintf("Calculating LD statistics for SNPs in %s pathway\n",
                    basename(file_path_sans_ext(bed[i])))) }

      } else {
        cat(sprintf("Reading input SNP set %s... ", basename(genoF)))
        start.time <- Sys.time()
        all.ss <- read.plink(genoF)
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)

        # Get n SNPs per high-confidence pathways - gives the number to
        # sample for each N random set (N = # of hc pathways)
        snps <- list.files(path=snpDir, pattern="*.snps$", full.names=T)
        hc.paths <- lapply(snps, readLines)
        hc.path.size <- sapply(hc.paths, length)

        set.seed(1)
        for (i in 1:length(hc.path.size)) {
          n <- hc.path.size[i]
          cat(sprintf("random set %i with %i snps\n", i, n))
          ss[[i]] <- all.ss$genotypes[, sample(ncol(all.ss$genotypes),
                                      n, replace=F)]
        }
    }
    for (i in 1:length(ss)) {

      if (pathSet != "rand") {
        ld.calc[[i]] <- ld(ss[[i]]$genotypes,
                           stats=c("D.prime", "R.squared"),
                           depth=ncol(ss[[i]]$genotypes)-1)
        snp.map <- ss[[i]]$map #genomic location of each SNP
      } else {
        ld.calc[[i]] <- ld(ss[[i]],
                           stats=c("D.prime", "R.squared"),
                           depth=ncol(ss[[i]])-1)
        snp.map <- all.ss$map }

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

  ### 1) Calc LD stats for high-confidence pathways SNPs ###
  cat("\nCalculating long-range correlation b/w high-confidence pathway SNPs.\n")

  calcLD(pathSet="high-confidence", snpDir=hcInDir)
  hc.all.pairs <- all.pairs
  hc.diff.num <- diff.num
  hc.diff.r2.mean <- diff.r2.mean
  hc.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  ### 2) Calc LD stats for low-confidence pathway SNPs ###
  Sys.sleep(5)
  cat("\nCalculating long-range correlation b/w low-confidence pathway SNPs.\n")

  calcLD(pathSet="low-confidence", snpDir=lcInDir)
  lc.all.pairs <- all.pairs
  lc.diff.num <- diff.num
  lc.diff.r2.mean <- diff.r2.mean
  lc.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  ### 3) Calc LD stats for random SNP sets ###
  Sys.sleep(5)
  cat("\nCalculating long-range correlation b/w randomly selected SNP sets.\n")

  calcLD(pathSet="rand", snpDir=hcInDir)
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
    hc <- data.frame(R.squared=hc.diff.r2.mean,
                     D.prime=hc.diff.dp.mean,
                     pathway.group="highconf")

    lc <- data.frame(R.squared=lc.diff.r2.mean,
                     D.prime=lc.diff.dp.mean,
                     pathway.group="lowconf")

    rn <- data.frame(R.squared=rn.diff.r2.mean,
                     D.prime=rn.diff.dp.mean,
                     pathway.group="rand")

    dist.dat <- rbind(hc, lc, rn)
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
        ggtitle(title) +
        theme_set(theme_minimal()) +
        theme(plot.title=element_text(hjust=0.5),
              text=element_text(size=17),
              legend.position="top",
              legend.title=element_blank(),
              panel.grid.major.x=element_blank())
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
        ggtitle(title) +
        theme_set(theme_minimal()) +
        theme(plot.title=element_text(hjust=0.5),
              text=element_text(size=17),
              legend.position="top",
              legend.title=element_blank(),
              panel.grid.major.x=element_blank())
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
          ggtitle(title) +
          theme_set(theme_minimal()) +
          theme(plot.title=element_text(hjust=0.5),
                text=element_text(size=17),
                legend.position="top",
                legend.title=element_blank(),
                panel.grid.major.x=element_blank(),
                axis.title.x=element_blank()) +
          scale_x_discrete(labels=paste("N=", table(dist.dat$pathway.group),
                                        sep=""))
      ggsave(sprintf("%s/boxplot_hc-lc-rn_r2.png", outDir),
          width=8, height=7.5)
  }
 else { return() }
}
