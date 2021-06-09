#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param inDir (char) path to files with pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @param repNum (interger) number of permutations to run
#' @param sampleNum (interger) number of SNPs to pick for each permutation
#' @return
#' @export

LDstats <- function(genoF, inDir, repNum=500L, sampleNum=400L,
                    makePlots=TRUE, outDir) {

  # Read PLINK files for high confidence pathway
  # NOTE: hc = high confidence
  hc.bed <- list.files(path=inDir, pattern="*.bed", full.names=T)
  hc.bim <- list.files(path=inDir, pattern="*.bim", full.names=T)
  hc.fam <- list.files(path=inDir, pattern="*.fam", full.names=T)

  hc.list <- c()
  hc.ld.calc <- list()
  hc.pairwise.df <- list()
  hc.diff.r2 <- list()
  hc.diff.dp <- list()

  message("\n-------HIGH CONFIDENCE PATHWAY SNPS-------\n")
  for (i in 1:length(hc.bed)) {

    # Convert PLINK files to snpStats input format
    # Output object is a list with 3 elements ($genotypes, $fam, $map)
    # NOTE: order is important!
    hc.list[[i]] <- read.plink(hc.bed[i], hc.bim[i], hc.fam[i])

    # Calculate linkage disequilbrium statistics
    # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
    # to be considered, so that depth=1 would specify calculation of LD b/w
    # immediately adjacent SNPs. For our purposes we want to determine LD
    # b/w all SNPs in each pathway despite their distance from each other,
    # so we specify depth as ((number of SNPs)-1)
    cat(sprintf("Calculating LD statistics for SNPs in %s pathway...",
                basename(file_path_sans_ext(hc.bed[i]))))

    hc.ld.calc[[i]] <- ld(hc.list[[i]]$genotypes,
                          stats=c("D.prime", "R.squared"),
		                      depth=ncol(hc.list[[i]]$genotypes)-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP-SNP pair
    snp.map <- hc.list[[i]]$map

    # Turn each LD matrix into a data frame
    hc.r2 <- as.matrix(hc.ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
    hc.r2 <- subset(melt(hc.r2), value!=0) #for all non-zero values
    colnames(hc.r2)[3] <- "R.squared"

    hc.dp <- as.matrix(hc.ld.calc[[i]]$D.prime)
    hc.dp <- subset(melt(hc.dp), value!=0)
    colnames(hc.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    hc.all.stats <- merge(hc.r2, hc.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(hc.all.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    hc.pairwise <- merge(snp.map, hc.all.stats, by="snp.name")
    colnames(hc.pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    hc.pairwise <- merge(snp.map, hc.pairwise, by="snp.name")
    colnames(hc.pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                               "chr_2", "pos_2", "R.squared", "D.prime")
    hc.pairwise$dist <- abs(hc.pairwise$pos_1 - hc.pairwise$pos_2)

    hc.pairwise.df[[i]] <- hc.pairwise %>% mutate(R.squared=round(R.squared, 3))
    hc.pairwise.df[[i]] <- hc.pairwise %>% mutate(D.prime=round(D.prime, 3))

    hc.diff.r2[[i]] <- filter(hc.pairwise.df[[i]], chr_1 != chr_2) %>%
                              dplyr::select(R.squared) %>% unlist
    hc.diff.dp[[i]] <- filter(hc.pairwise.df[[i]], chr_1 != chr_2) %>%
                              dplyr::select(D.prime) %>% unlist

    cat(" done.\n")
  }

  hc.all <- do.call("rbind", hc.pairwise.df)
  hc.diff.num <- sapply(hc.diff.r2, length) #get sample size per pathway
  hc.diff.r2.mean <- sapply(hc.diff.r2, mean)   #mean r2 value per pathway
  hc.diff.dp.mean <- sapply(hc.diff.dp, mean) #mean dprime per pathway

  cat(sprintf("\nCalculated LD for %i total SNP-SNP pairs.\n", nrow(hc.all)))
  cat(sprintf("\t --> %i total interchromosomal pairs.\n", sum(hc.diff.num)))

  #remove original data objects to clear memory
  rm(hc.list, hc.ld.calc)

  #============================================================================#
  # Permute random samples from original PLINK genotype data and calculate LD
  message("\n-------RANDOMLY SELECTED SNPS-------\n")

  # Large vector, time intensive depending on size of file
  cat(sprintf("Reading input SNP set %s... ", basename(genoF)))
  start.time <- Sys.time()
  test <- read.plink(genoF)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken) #print time taken to read in PLINK files

  null.ld.calc <- list()
  null.pairwise.df <- list()
  null.diff.r2 <- list()
  null.diff.dp <- list()

  Sys.sleep(5)
  for (i in 1:repNum) {
    cat(sprintf("Calculating LD within random sample matrix %i...", i))

    # Generating LD stats for 500 permutations of 400 SNPs each
    # later will plot mean null r2 / dprime distribution per perm via ggplot
    null.ld.calc[[i]] <- ld(test$genotypes[, sample(ncol(test$genotypes),
                                          sampleNum, replace=F)],
                            stats=c("D.prime", "R.squared"),
                            depth=sampleNum-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP pair
    snp.map <- test$map

    # Turn each LD matrix into a data frame
    null.r2 <- as.matrix(null.ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
    null.r2 <- subset(melt(null.r2), value!=0) #melt df and remove '0's
    colnames(null.r2)[3] <- "R.squared"

    null.dp <- as.matrix(null.ld.calc[[i]]$D.prime)
    null.dp <- subset(melt(null.dp), value!=0)
    colnames(null.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    null.stats <- merge(null.r2, null.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(null.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    null.pairwise <- merge(snp.map, null.stats, by="snp.name")
    colnames(null.pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    null.pairwise <- merge(snp.map, null.pairwise, by="snp.name")
    colnames(null.pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                                 "chr_2", "pos_2", "R.squared", "D.prime")
    # Calculate distance between SNP pairs
    null.pairwise$dist <- abs(null.pairwise$pos_1 - null.pairwise$pos_2)

    # Round r2 and Dprime values to 3 decimal points
    null.pairwise.df[[i]] <- null.pairwise %>% mutate(R.squared=round(R.squared, 3))
    null.pairwise.df[[i]] <- null.pairwise %>% mutate(D.prime=round(D.prime, 3))

    # Used to build null distruibution of mean R2 / Dprime for SNP-SNP pairs on
    # different chromosomes per sample 'pathway'
    null.diff.r2[[i]] <- filter(null.pairwise.df[[i]], chr_1 != chr_2) %>%
                                dplyr::select(R.squared) %>% unlist
    null.diff.dp[[i]] <- filter(null.pairwise.df[[i]], chr_1 != chr_2) %>%
                                dplyr::select(D.prime) %>% unlist

    cat(" done.\n")
  }

#  null.all <- do.call("rbind", null.pairwise.df)
  null.diff.num <- sapply(null.diff.r2, length)
  null.diff.r2.mean <- sapply(null.diff.r2, mean)
  null.diff.dp.mean <- sapply(null.diff.dp, mean)

#  cat(sprintf("Calculated LD for %i randomly permuted SNP interactions.\n",
#              nrow(null.all)))

  rm(test, null.ld.calc)

  #============================================================================#
  ## PLOT STATS
  if (makePlots == TRUE) {
    message("\nGenerating LD statistics plots.\n")

    # Set variables and other functions
    title <- paste("Degree of co-selection per interchromosomal SNP-SNP",
                   "\ninteraction within the high-confidence pathways",
                   "\nvs. randomly selected SNP groups")
    r2.axis.title <- bquote("Pairwise LD value (mean " *r^2*" per pathway)")
    dp.axis.title <- bquote("Pairwise LD value (mean D' per pathway)")

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
    calc.p <- function(null, real){
      return( (length(which(null > mean(real)))+1) / (length(null)+1) )
    }

    #################### PLOT 1: Null distribution vs real #####################
    # Plotting mean R2 value for random SNP-SNP pairs on diff chromosome
    # against real mean R2 value for high-conf SNPs
    # interchromosomal SNP pairs used as proxy for co-selection/genetic ixns
    real <- data.frame(R.squared=hc.diff.r2.mean,
                       D.prime=hc.diff.dp.mean,
                       pathway.group="highconf")

    null <- data.frame(R.squared=null.diff.r2.mean,
                       D.prime=null.diff.dp.mean,
                       pathway.group="rand")

    dist.dat <- rbind(real, null)
    mean.dat <- ddply(dist.dat, "pathway.group", summarise,
                      R.squared.mean=mean(R.squared),
                      D.prime.mean=mean(D.prime))

    ggplot(dist.dat, aes(x=R.squared, colour=pathway.group,
                         fill=pathway.group)) +
        #geom_density(alpha=0.3) +
        geom_histogram(bins=20, alpha=0.5, position="identity") +
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
    ggsave(sprintf("%s/dist_nullvreal_r2.png", outDir),
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
    ggsave(sprintf("%s/dist_nullvreal_r2_facet.png", outDir),
        width=8, height=7)

    perm.p = calc.p(null.diff.r2.mean, hc.diff.r2.mean)

    cat("\n-----------------------------------------------------------------")
    cat(sprintf("\nP-value for permuted sample vs. real test statistic = %g\n",
                perm.p))
       if (perm.p < 0.05) {
         cat("\nCelebrate! It's significant.")
       } else { cat("\nThat doesn't look so good.") }
    cat("\n-----------------------------------------------------------------\n")

    ############################ PLOT 2: Boxplots ##############################
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
      ggsave(sprintf("%s/boxplot_nullvsreal_r2.png", outDir),
          width=8, height=7.5)

    ####################### PLOT 3: Boxplot per pathway ########################
    # Boxplot of LD stats per high confidence pathway
    hc.r2.df <- melt(hc.diff.r2, value.name="R.squared")
    hc.dp.df <- melt(hc.diff.dp, value.name="D.prime")

    title <- paste("Degree of co-selection per interchromosomal SNP-SNP",
                   "\npair within each high-confidence pathway")
    axis.title2 <- "# of interchromosomal SNP-SNP pairs per pathway"

    ggplot(hc.r2.df, aes(x=factor(L1), y=R.squared,
                         color="#F8766D", fill="#F8766D")) +
         geom_boxplot() +
         stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
                      fun.data=box.style) +
         stat_summary(geom="text", color="white", fun.data=give.n,
                      position=position_dodge(width=0.75)) +
         scale_y_continuous(r2.axis.title) +
         scale_x_discrete(axis.title2)+
         ggtitle(title) +
         theme(axis.text.x=element_text(vjust=0.4, hjust=1)) +
         theme_set(theme_minimal()) +
         theme(plot.title=element_text(hjust=0.5),
               text=element_text(size=17),
               legend.position="none",
               panel.grid.major.x=element_blank()) +
         geom_hline(yintercept=0.1, colour="grey", linetype="dashed", size=1)
    ggsave(sprintf("%s/hc_bars_r2.png", outDir), width=11)
  }
 else { return() }
}
