#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param lcInDir (char) path to files with low-confidence pathway SNP lists
#' @param psInDir (char) path to files with pseudo pathway SNP lists
#' @param allInDir (char) path to files with all pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @param outDir (char) path to output directory
#' @return
#' @export

LDstats <- function(hcInDir, lcInDir, psInDir, allInDir,
                    makePlots=TRUE, outDir) {

  # Calculate linkage disequilbrium statistics
  # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
  # to be considered, so that depth=1 would specify calculation of LD b/w
  # immediately adjacent SNPs. For our purposes we want to determine LD
  # b/w all SNPs in each pathway despite their distance from each other,
  # so we specify depth in ld() as ((number of SNPs)-1)
  calcLD <- function(snpDir, popcode) {

  ss <- c(); snp.map <- c()
  ld.calc <- list(); pairwise.df <- list(); diff <- list()
  diff.r2 <- list(); diff.dp <- list()

  # Read PLINK files for each pathway
  bed <- list.files(path=snpDir, pattern="*.bed", full.names=T)
  bim <- list.files(path=snpDir, pattern="*.bim", full.names=T)
  fam <- list.files(path=snpDir, pattern="*.fam", full.names=T)

  # Convert PLINK files to snpStats input format
  # Output object is a list with 3 elements ($genotypes, $fam, $map)
  # NOTE: order is important!
  for (i in 1:length(bed)) {
    cat(sprintf("*Reading pathway %s PLINK set\n",
          basename(file_path_sans_ext(bed[i]))))

    ss[[i]] <- read.plink(bed[i], bim[i], fam[i])

    ############################
    # experimenting with pop-based ld calc
    if (popcode == 0) {
      pop <- which(ss[[i]]$fam$affected != popcode)
      ss[[i]]$genotypes <- ss[[i]]$genotypes[pop, ]

      cat(sprintf("*Keeping %i rows in SnpMatrix for both population samples\n",
              length(pop)))
      print(ss[[i]]$genotypes)
    } else {
      pop <- which(ss[[i]]$fam$affected == popcode)
      ss[[i]]$genotypes <- ss[[i]]$genotypes[pop, ]

      cat(sprintf("*Subsetting %i genotypes based on %s population(s)\n",
          length(pop),
          if (popcode == 1) {pop1}
          else if (popcode == 2) {pop2} ))
      print(ss[[i]]$genotypes)
    }
    ############################

    cat("*Calculating LD statistics\n")
    ld.calc[[i]] <- ld(ss[[i]]$genotypes,
                       stats=c("D.prime", "R.squared"),
                       depth=ncol(ss[[i]]$genotypes)-1)
    snp.map <- ss[[i]]$map #genomic location of each SNP

    # Turn each LD matrix into a data frame
    r2 <- as.matrix(ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
    r2 <- subset(melt(r2), value!=0) #keep all non-zero values
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

    cat("*Subsetting for interchromosomal SNP-SNP pairs\n")
    diff[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2)

    diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
    diff.dp[[i]] <- select(diff[[i]], D.prime) %>% unlist
    cat("done.\n\n")
  }

  all.pairs <<- do.call("rbind", pairwise.df)
  diff.pairs <<- do.call("rbind", diff)
  diff.num <<- sapply(diff.r2, length)   #all SNP-SNP pairs per path
  diff.r2.mean <<- sapply(diff.r2, mean) #mean r2 value per pathway
  diff.dp.mean <<- sapply(diff.dp, mean) #mean dprime per pathway

  cat(sprintf("\nCalculated LD for %i total SNP pairs.\n", nrow(all.pairs)))
  cat(sprintf("\t --> %i total interchromosomal pairs.\n", sum(diff.num)))
}

  # How all pairwise combinations are calculated
  # n is the number of SNPs in each pathway, divided by 2 to account for
  # duplicate pairs ("n choose 2")
  combos <- function(n) {
    return( n*(n-1)/2 ) }

  ### 1) Calc LD stats for high-confidence pathways SNPs ###
  cat("\nCalculating long-range correlation b/w high-confidence pathway SNPs.\n")

  calcLD(snpDir=hcInDir, popcode=0)
  hc.all.pairs <- all.pairs
  hc.diff.pairs <<- diff.pairs
  hc.diff.num <- diff.num
  hc.diff.r2.mean <- diff.r2.mean
  hc.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  ### 2) Calc LD stats for low-confidence pathway SNPs ###
  Sys.sleep(5)
  cat("\nCalculating long-range correlation b/w low-confidence pathway SNPs.\n")

  calcLD(snpDir=lcInDir, popcode=0)
  lc.all.pairs <- all.pairs
  lc.diff.pairs <<- diff.pairs
  lc.diff.num <- diff.num
  lc.diff.r2.mean <- diff.r2.mean
  lc.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  ### 3) Calc LD stats for random SNP sets ###
  Sys.sleep(5)
  cat("\nCalculating long-range correlation b/w randomly selected SNP sets.\n")

  calcLD(snpDir=pseudoInDir)
  ps.all.pairs <- all.pairs
  ps.diff.num <- diff.num
  ps.diff.r2.mean <- diff.r2.mean
  ps.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

  calcLD(snpDir=allInDir)
  all.all.pairs <- all.pairs
  all.diff.num <- diff.num
  all.diff.r2.mean <- diff.r2.mean
  all.diff.dp.mean <- diff.dp.mean
  rm(all.pairs, diff.num, diff.r2.mean, diff.dp.mean)

### DEBUGGING
all.hc.snps <- c(hc.all.pairs$snp_1, hc.all.pairs$snp_2)
unique.hc.snps <- as.data.frame(unique(all.hc.snps))
colnames(unique.hc.snps) <- "snp"

all.lc.snps <- c(lc.all.pairs$snp_1, lc.all.pairs$snp_2)
unique.lc.snps <- as.data.frame(unique(all.lc.snps))
colnames(unique.lc.snps) <- "snp"

all.ps.snps <- c(ps.all.pairs$snp_1, ps.all.pairs$snp_2)
unique.ps.snps <- as.data.frame(unique(all.ps.snps))
colnames(unique.ps.snps) <- "snp"

hc.lc <- merge(unique.hc.snps, unique.lc.snps, by="snp")
hc.ps <- merge(unique.hc.snps, unique.ps.snps, by="snp")
lc.ps <- merge(unique.lc.snps, unique.ps.snps, by="snp")

nrow(hc.lc)
#[1] 106

nrow(hc.ps)
#[1] 487  too many hc snps in pseudo paths, control for?

nrow(lc.ps)
#[1] 226

  #============================================================================#
  ## PLOT STATS
  if (makePlots == TRUE) {
    message("\nGenerating LD statistics plots.\n")

    # Set variables and other functions
    title <- "Degree of co-selection per inter-chromosomal\nSNP-SNP interaction"
    r2.axis.title <- bquote(bold("Pairwise LD value per SNP-SNP interaction"))
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
    calc.p <- function(rn, hc){
      return( (length(which(ps > mean(hc)))+1) / (length(ps)+1) )
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

    #hc.ceu <- data.frame(R.squared=hc.diff.r2.mean.ceu,
    #                     D.prime=hc.diff.dp.mean.ceu,
    #                     pathway.group="highconf",
    #                     population="European")

    #lc.ceu <- data.frame(R.squared=lc.diff.r2.mean.ceu,
    #                     D.prime=lc.diff.dp.mean.ceu,
    #                     pathway.group="lowconf",
    #                     population="European")

    #hc.asw <- data.frame(R.squared=hc.diff.r2.mean.asw,
    #                     D.prime=hc.diff.dp.mean.asw,
    #                     pathway.group="highconf",
    #                     population="African")

    #lc.asw <- data.frame(R.squared=lc.diff.r2.mean.asw,
    #                     D.prime=lc.diff.dp.mean.asw,
    #                     pathway.group="lowconf",
    #                     population="African")

  ######## w/all r2 values (not just means)
  hc.diff.pairs$pathway.group <- "Enriched"
  lc.diff.pairs$pathway.group <- "Nonenriched"
  dat <- rbind(hc.diff.pairs, lc.diff.pairs)

  mean.dat <- ddply(dat, "pathway.group", summarise,
                    R.squared.mean=mean(R.squared))

  ### sig via KS test
  pval <- ks.test(hc.diff.pairs$R.squared, lc.diff.pairs$R.squared,
                  alternative="less")
  p_text <- pval$p.value

  # density distribution plot
  p1 <- ggplot(dat, aes(x=R.squared, colour=pathway.group,
               fill=pathway.group)) +
          #  facet_grid(population ~ .) +
            geom_density(alpha=0.2) +
          #  geom_vline(data=mean.dat, aes(xintercept=R.squared.mean,
          #                                colour=pathway.group),
          #             linetype="dashed", size=1) +
            scale_x_continuous("") +
            scale_y_continuous("Density") +
            ggtitle(title) +
          #  xlim(0.2,1) +
            theme_Publication() +
            scale_colour_Publication() +
            scale_fill_Publication() +
            annotate("text", x=0.07, y=25,
                     label=sprintf("paste(italic(P), \" = %1.3g\")", p_text),
                     parse=TRUE)

    # eCDF plot (cumulative density at each r2)
    p2 <- ggplot(dat, aes(x=R.squared, colour=pathway.group)) +
            stat_ecdf() +
            scale_x_continuous(r2.axis.title) +
            scale_y_continuous("Cumulative density") +
            theme_Publication() +
            scale_colour_Publication() +
            scale_fill_Publication()

    both <- arrangeGrob(p1, p2, ncol=1)
    ggsave("density_ecdf_hc_lc-nes_r2.png", both, width=6.5, height=8)

    # Boxplot of total LD stat distribution b/w null vs. real
#p<-  ggplot(dist.dat, aes(x=pathway.group, y=R.squared)) +
#          facet_grid(population ~ .) +
#          geom_boxplot(outlier.colour=NULL,
#                       aes(colour=pathway.group, fill=pathway.group)) +
#          stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
#                       fun.data=box.style) +
#          scale_y_continuous(r2.axis.title) +
#          ggtitle(title) +
#          theme_set(theme_minimal()) +
#          theme(plot.title=element_text(hjust=0.5),
#                text=element_text(size=17),
#                legend.position="top",
#                legend.title=element_blank(),
#                panel.grid.major.x=element_blank(),
#                axis.title.x=element_blank()) +
#          scale_x_discrete(labels=paste("N=", table(dist.dat$pathway.group),
#                                        sep=""))
#      ggsave(sprintf("%s/boxplot_hc-lc-ps-rn-all_r2_pop-strat.png", outDir),p,
#          width=8, height=7.5)
#  }
 else { return() }
}
