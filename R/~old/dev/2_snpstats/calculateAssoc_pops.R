#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param plinkF (char) path to file with SNP genotype data (PLINK format)
#' @param highSnpDir (char) path to files with pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @return
#' @export

calculateAssoc_pops <- function(plinkPop1F, plinkPop2,
                                hcPop1Dir, hcPop2Dir,
                                makePlots=FALSE) {

  # Read PLINK files for each high confidence pathway
  # NOTE: hc = high confidence
  hc.bed.pop1 <- list.files(path=hcPop1Dir, pattern="*.bed", full.names=T)
  hc.bim.pop1 <- list.files(path=hcPop1Dir, pattern="*.bim", full.names=T)
  hc.fam.pop1 <- list.files(path=hcPop1Dir, pattern="*.fam", full.names=T)

  hc.list.pop1 <- c()
  hc.ld.calc.pop1 <- list()
  hc.pair.df.pop1 <- list()
  hc.diff.r2.pop1 <- list()
  hc.diff.dp.pop1 <- list()

  message("\n-------HIGH CONFIDENCE PATHWAY SNPS -- POPULATION 1 -------\n")
  for (i in 1:length(hc.bed.pop1)) {

    # Convert PLINK files to snpStats input format
    # Output object is a list with 3 elements ($genotypes, $fam, $map)
    # NOTE: order is important!
    hc.list.pop1[[i]] <- read.plink(hc.bed.pop1[i], hc.bim.pop1[i],
                                    hc.fam.pop1[i])

    # Calculate linkage disequilbrium statistics
    # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
    # to be considered, so that depth=1 would specify calculation of LD b/w
    # immediately adjacent SNPs. For our purposes we want to determine LD
    # b/w all SNPs in each pathway despite their distance from each other,
    # so we specify depth as ((number of SNPs)-1)
    cat(sprintf("Calculating LD statistics for SNPs in %s pathway...",
                basename(file_path_sans_ext(hc.bed.pop1[i]))))

    hc.ld.calc.pop1[[i]] <- ld(hc.list.pop1[[i]]$genotypes,
                               stats=c("D.prime", "R.squared"),
		                           depth=ncol(hc.list.pop1[[i]]$genotypes)-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP-SNP pair
    snp.map <- hc.list.pop1[[i]]$map

    # Turn each LD matrix into a data frame
    hc.r2 <- as.matrix(hc.ld.calc.pop1[[i]]$R.squared) #convert sparseMatrix to regular matrix
    hc.r2 <- subset(melt(hc.r2), value!=0) #for all non-zero values
    colnames(hc.r2)[3] <- "R.squared"

    hc.dp <- as.matrix(hc.ld.calc.pop1[[i]]$D.prime)
    hc.dp <- subset(melt(hc.dp), value!=0)
    colnames(hc.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    hc.all.stats <- merge(hc.r2, hc.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(hc.all.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    hc.pair <- merge(snp.map, hc.all.stats, by="snp.name")
    colnames(hc.pair)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    hc.pair <- merge(snp.map, hc.pair, by="snp.name")
    colnames(hc.pair) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                           "chr_2", "pos_2", "R.squared", "D.prime")
    hc.pair$dist <- abs(hc.pair$pos_1 - hc.pair$pos_2)

    hc.pair.df.pop1[[i]] <- hc.pair %>% mutate(R.squared=round(R.squared, 3))
    hc.pair.df.pop1[[i]] <- hc.pair %>% mutate(D.prime=round(D.prime, 3))

    hc.diff.r2.pop1[[i]] <- filter(hc.pair.df.pop1[[i]], chr_1 != chr_2) %>%
                            select(R.squared) %>% unlist
    hc.diff.dp.pop1[[i]] <- filter(hc.pair.df.pop1[[i]], chr_1 != chr_2) %>%
                            select(D.prime) %>% unlist

    cat(" done.\n")
  }

  # Same calculations as above for population 2
  hc.bed.pop2 <- list.files(path=hcPop2Dir, pattern="*.bed", full.names=T)
  hc.bim.pop2 <- list.files(path=hcPop2Dir, pattern="*.bim", full.names=T)
  hc.fam.pop2 <- list.files(path=hcPop2Dir, pattern="*.fam", full.names=T)

  hc.list.pop2 <- c()
  hc.ld.calc.pop2 <- list()
  hc.pair.df.pop2 <- list()
  hc.diff.r2.pop2 <- list()
  hc.diff.dp.pop2 <- list()

  message("\n-------HIGH CONFIDENCE PATHWAY SNPS -- POPULATION 2 -------\n")
  for (i in 1:length(hc.bed.pop2)) {

    # Convert PLINK files to snpStats input format
    # Output object is a list with 3 elements ($genotypes, $fam, $map)
    # NOTE: order is important!
    hc.list.pop2[[i]] <- read.plink(hc.bed.pop2[i], hc.bim.pop2[i],
                                    hc.fam.pop2[i])

    # Calculate linkage disequilbrium statistics (R squared)
    # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
    # to be considered, so that depth=1 would specify calculation of LD b/w
    # immediately adjacent SNPs. For our purposes we want to determine LD
    # b/w all SNPs in each pathway despite their distance from each other,
    # so we specify depth as ((number of SNPs)-1)
    cat(sprintf("Calculating LD statistics for SNPs in %s pathway...",
                basename(file_path_sans_ext(hc.bed.pop2[i]))))

    hc.ld.calc.pop2[[i]] <- ld(hc.list.pop2[[i]]$genotypes,
                               stats=c("D.prime", "R.squared"),
		                           depth=ncol(hc.list.pop2[[i]]$genotypes)-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP-SNP pair
    snp.map <- hc.list.pop2[[i]]$map

    # Turn each LD matrix into a data frame
    hc.r2 <- as.matrix(hc.ld.calc.pop2[[i]]$R.squared) #convert sparseMatrix to regular matrix
    hc.r2 <- subset(melt(hc.r2), value!=0) #for all non-zero values
    colnames(hc.r2)[3] <- "R.squared"

    hc.dp <- as.matrix(hc.ld.calc.pop2[[i]]$D.prime)
    hc.dp <- subset(melt(hc.dp), value!=0)
    colnames(hc.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    hc.all.stats <- merge(hc.r2, hc.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(hc.all.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    hc.pair <- merge(snp.map, hc.all.stats, by="snp.name")
    colnames(hc.pair)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    hc.pair <- merge(snp.map, hc.pair, by="snp.name")
    colnames(hc.pair) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                           "chr_2", "pos_2", "R.squared", "D.prime")
    hc.pair$dist <- abs(hc.pair$pos_1 - hc.pair$pos_2)

    hc.pair.df.pop2[[i]] <- hc.pair %>% mutate(R.squared=round(R.squared, 3))
    hc.pair.df.pop2[[i]] <- hc.pair %>% mutate(D.prime=round(D.prime, 3))

    hc.diff.r2.pop2[[i]] <- filter(hc.pair.df.pop2[[i]], chr_1 != chr_2) %>%
                            select(R.squared) %>% unlist
    hc.diff.dp.pop2[[i]] <- filter(hc.pair.df.pop2[[i]], chr_1 != chr_2) %>%
                            select(D.prime) %>% unlist

    cat(" done.\n")
  }

  hc.all.pop1 <- do.call("rbind", hc.pair.df.pop1)
  hc.diff.num.pop1 <- sapply(hc.diff.r2.pop1, length) #sample size per pathway
  hc.diff.r2.mean.pop1 <- sapply(hc.diff.r2.pop1, mean) #mean r2 per pathway
  hc.diff.dp.mean.pop1 <- sapply(hc.diff.dp.pop1, mean) #mean dprime per pathway

  hc.all.pop2 <- do.call("rbind", hc.pair.df.pop2)
  hc.diff.num.pop2 <- sapply(hc.diff.r2.pop2, length)
  hc.diff.r2.mean.pop2 <- sapply(hc.diff.r2.pop2, mean)
  hc.diff.dp.mean.pop2 <- sapply(hc.diff.dp.pop2, mean)

  #remove original data objects to clear memory
  rm(hc.list.pop1, hc.list.pop2, hc.ld.calc.pop1, hc.ld.calc.pop2)

  #============================================================================#
  # Permute random samples from original PLINK genotype data and calculate LD
  message("\n-------RANDOMLY SELECTED SNPS-------\n")

  rep.num <- 500L      #how many permutations to run
  sample.num <- 400L   #number of SNPs to pick for each permutation

  # Large vector, time intensive depending on size of file
  start.time <- Sys.time()
  pop1 <- read.plink(plinkPop1F)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken) #print time taken to read in PLINK files

  null.ld.calc.pop1 <- list()
  null.pair.df.pop1 <- list()
  null.diff.r2.pop1 <- list()
  null.diff.dp.pop1 <- list()

  for (i in 1:rep.num) {
    cat(sprintf("Calculating LD within random sample matrix %i...", i))

    # Generating LD stats for 500 permutations of 400 SNPs each
    # later will plot mean null r2 / dprime distribution per perm via ggplot
    null.ld.calc.pop1[[i]] <- ld(pop1$genotypes[, sample(ncol(pop1$genotypes),
                                 sample.num, replace=F)],
                                 stats=c("D.prime", "R.squared"),
                                 depth=sample.num-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP pair
    snp.map <- pop1$map

    # Turn each LD matrix into a data frame
    null.r2 <- as.matrix(null.ld.calc.pop1[[i]]$R.squared) #convert sparseMatrix to regular matrix
    null.r2 <- subset(melt(null.r2), value!=0) #melt df and remove '0's
    colnames(null.r2)[3] <- "R.squared"

    null.dp <- as.matrix(null.ld.calc.pop1[[i]]$D.prime)
    null.dp <- subset(melt(null.dp), value!=0)
    colnames(null.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    null.stats <- merge(null.r2, null.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(null.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    null.pair <- merge(snp.map, null.stats, by="snp.name")
    colnames(null.pair)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    null.pair <- merge(snp.map, null.pair, by="snp.name")
    colnames(null.pair) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                             "chr_2", "pos_2", "R.squared", "D.prime")
    # Calculate distance between SNP pairs
    null.pair$dist <- abs(null.pair$pos_1 - null.pair$pos_2)

    # Round r2 and Dprime values to 3 decimal points
    null.pair.df.pop1[[i]] <- null.pair %>% mutate(R.squared=round(R.squared, 3))
    null.pair.df.pop1[[i]] <- null.pair %>% mutate(D.prime=round(D.prime, 3))

    # Used to build null distruibution of mean R2 / Dprime for SNP-SNP pairs on
    # different chromosomes per sample 'pathway'
    null.diff.r2.pop1[[i]] <- filter(null.pair.df.pop1[[i]], chr_1 != chr_2) %>%
                              select(R.squared) %>% unlist
    null.diff.dp.pop1[[i]] <- filter(null.pair.df.pop1[[i]], chr_1 != chr_2) %>%
                              select(D.prime) %>% unlist

    cat(" done.\n")
  }

  start.time <- Sys.time()
  pop2 <- read.plink(plinkPop2F)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken) #print time taken to read in PLINK files

  null.ld.calc.pop2 <- list()
  null.pair.df.pop2 <- list()
  null.diff.r2.pop2 <- list()
  null.diff.dp.pop2 <- list()

  for (i in 1:rep.num) {
    cat(sprintf("Calculating LD within random sample matrix %i...", i))

    # Generating LD stats for 500 permutations of 400 SNPs each
    # later will plot mean null r2 / dprime distribution per perm via ggplot
    null.ld.calc.pop2[[i]] <- ld(pop2$genotypes[, sample(ncol(pop2$genotypes),
                                 sample.num, replace=F)],
                                 stats=c("D.prime", "R.squared"),
                                 depth=sample.num-1)

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP pair
    snp.map <- pop2$map

    # Turn each LD matrix into a data frame
    null.r2 <- as.matrix(null.ld.calc.pop2[[i]]$R.squared) #convert sparseMatrix to regular matrix
    null.r2 <- subset(melt(null.r2), value!=0) #melt df and remove '0's
    colnames(null.r2)[3] <- "R.squared"

    null.dp <- as.matrix(null.ld.calc.pop2[[i]]$D.prime)
    null.dp <- subset(melt(null.dp), value!=0)
    colnames(null.dp)[3] <- "D.prime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    null.stats <- merge(null.r2, null.dp, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(null.stats)[1] <- "snp.name"
    snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

    null.pair <- merge(snp.map, null.stats, by="snp.name")
    colnames(null.pair)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    null.pair <- merge(snp.map, null.pair, by="snp.name")
    colnames(null.pair) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                             "chr_2", "pos_2", "R.squared", "D.prime")
    # Calculate distance between SNP pairs
    null.pair$dist <- abs(null.pair$pos_1 - null.pair$pos_2)

    # Round r2 and Dprime values to 3 decimal points
    null.pair.df.pop2[[i]] <- null.pair %>% mutate(R.squared=round(R.squared, 3))
    null.pair.df.pop2[[i]] <- null.pair %>% mutate(D.prime=round(D.prime, 3))

    # Used to build null distruibution of mean R2 / Dprime for SNP-SNP pairs on
    # different chromosomes per sample 'pathway'
    null.diff.r2.pop2[[i]] <- filter(null.pair.df.pop2[[i]], chr_1 != chr_2) %>%
                              select(R.squared) %>% unlist
    null.diff.dp.pop2[[i]] <- filter(null.pair.df.pop2[[i]], chr_1 != chr_2) %>%
                              select(D.prime) %>% unlist

    cat(" done.\n")
  }

  #null.all <- do.call("rbind", null.pairwise.df)
  null.diff.num.pop1 <- sapply(null.diff.r2.pop1, length)
  null.diff.r2.mean.pop1 <- sapply(null.diff.r2.pop1, mean)
  null.diff.dp.mean.pop1 <- sapply(null.diff.dp.pop1, mean)

  null.diff.num.pop2 <- sapply(null.diff.r2.pop2, length)
  null.diff.r2.mean.pop2 <- sapply(null.diff.r2.pop2, mean)
  null.diff.dp.mean.pop2 <- sapply(null.diff.dp.pop2, mean)

  rm(pop1, pop2, null.ld.calc.pop1, null.ld.calc.pop2)

  #============================================================================#
  ## PLOT STATS
  message("\n-------PLOTS-------\n")
  cat("Linkage disequilbrium statistic plots...")

  # Set variables and other functions
  title <- paste("Degree of co-selection per interchromosomal SNP-SNP",
                 "\ninteraction within the high-confidence pathways",
                 "\nvs. randomly selected SNP groups")
  r2.axis.title <- bquote("Pairwise LD value (mean " *r^2*" per pathway)")
  dp.axis.title <- bquote("Pairwise LD value (mean D' per pathway)")

  ## for use with stat_summary(fun.data=box.style); allows white median line to
  ## appear after colouring and filling boxplots
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

  #################### PLOT 1: Null distribution vs real #######################
  # Plotting mean R2 value for random SNP-SNP pairs on diff chromosome
  # against real mean R2 value for high-conf SNPs
  # interchromosomal SNP pairs used as proxy for co-selection/genetic ixns
  real.pop1 <- data.frame(R.squared=hc.diff.r2.mean.pop1,
                          D.prime=hc.diff.dp.mean.pop1,
                          pathway.group="highconf")

  null.pop1 <- data.frame(R.squared=null.diff.r2.mean.pop1,
                          D.prime=null.diff.dp.mean.pop1,
                          pathway.group="rand")

  dist.dat.pop1 <- rbind(real.pop1, null.pop1)
  mean.dat.pop1 <- ddply(dist.dat.pop1, "pathway.group", summarise,
                         R.squared.mean=mean(R.squared),
                         D.prime.mean=mean(D.prime))
 ####
 real.pop2 <- data.frame(R.squared=hc.diff.r2.mean.pop2,
                         D.prime=hc.diff.dp.mean.pop2,
                         pathway.group="highconf")

 null.pop2 <- data.frame(R.squared=null.diff.r2.mean.pop2,
                         D.prime=null.diff.dp.mean.pop2,
                         pathway.group="rand")

 dist.dat.pop2 <- rbind(real.pop2, null.pop1)
 mean.dat.pop2 <- ddply(dist.dat.pop2, "pathway.group", summarise,
                        R.squared.mean=mean(R.squared),
                        D.prime.mean=mean(D.prime))


  ggplot(dist.dat.pop1, aes(x=R.squared, colour=pathway.group,
                            fill=pathway.group)) +
      geom_density(alpha=0.3) +
      #geom_histogram(bins=20, alpha=0.5, position="identity") +
      geom_vline(data=mean.dat.pop1, aes(xintercept=R.squared.mean,
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
  ggsave("dist_nullvreal_r2_ceu.png", width=8, height=7)

  perm.p = calc.p(null.diff.r2.mean, hc.diff.r2.mean)
  cat(sprintf("P-value for permuted sample vs. real test statistic= %g\n",
              perm.p))

  ############################ PLOT 2: Boxplots ################################
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
    ggsave("boxplot_nullvsreal_r2.png", width=8, height=7.5)

  ####################### PLOT 3: Boxplot per pathway ##########################
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
  ggsave("hc_bars_r2.png", width=11)

}
