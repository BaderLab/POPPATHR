#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param plinkF (char) path to file with SNP genotype data (PLINK format)
#' @param highSnpDir (char) path to files with pathway SNP lists
#' @param makePlots (logical) set to TRUE to generate plots
#' @return
#' @export

calculateAssoc <- function(plinkF, highSnpDir, makePlots=FALSE) {

  # Read PLINK files for high confidence pathway
  top_bed <- list.files(path=highSnpDir, pattern="*.bed", full.names=T)
  top_bim <- list.files(path=highSnpDir, pattern="*.bim", full.names=T)
  top_fam <- list.files(path=highSnpDir, pattern="*.fam", full.names=T)

  top_list <- c()
  top_r2_matrix <- list()
  top_dprime_matrix <- list()
  top_pairwise_df <- list()
  top_diff_chr <- list()

  message("\n-------HIGH CONFIDENCE PATHWAY SNPS-------\n")
  for (i in 1:length(top_bed)) {

    # Convert PLINK files to snpStats input format
    # Output object is a list with 3 elements ($genotypes, $fam, $map)
    # NOTE: order is important!
    top_list[[i]] <- read.plink(top_bed[i], top_bim[i], top_fam[i])

    # Calculate linkage disequilbrium statistics (R squared)
    # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
    # to be considered, so that depth=1 would specify calculation of LD b/w
    # immediately adjacent SNPs. For our purposes we want to determine LD
    # b/w all SNPs in each pathway despite their distance from each other,
    # so we specify depth as ((number of SNPs)-1)
    cat(sprintf("Calculating LD statistics for SNPs in %s pathway...",
                basename(file_path_sans_ext(top_bed[i]))))

    top_ld_calc[[i]] <- ld(top_list[[i]]$genotypes,
                           stats=c("D.prime", "R.squared"),
		                       depth=ncol(top_list[[i]]$genotypes)-1)
    top_r2_matrix[[i]] <- top_ld_calc[[i]]$R.squared
    top_dprime_matrix[[i]] <- top_ld_calc[[i]]$D.prime
    cat(" done.\n")

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP pair
    snp_map <- top_list[[i]]$map

    # Turn each LD matrix into a data frame
    top_r2 <- as.matrix(top_r2_matrix[[i]]) #convert sparseMatrix to regular matrix
    top_r2 <- subset(melt(top_r2), value!=0) #for all non-zero values
    colnames(top_r2)[3] <- "R2"

    top_dprime <- as.matrix(top_dprime_matrix[[i]])
    top_dprime <- subset(melt(top_dprime), value!=0)
    colnames(top_dprime)[3] <- "Dprime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    top_stats <- merge(top_r2, top_dprime, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(top_stats)[1] <- "snp.name"
    snp_map <- subset(snp_map, select=c("snp.name", "chromosome", "position"))
    top_pairwise <- merge(snp_map, top_stats, by="snp.name")
    colnames(top_pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    top_pairwise <- merge(snp_map, top_pairwise, by="snp.name")
    colnames(top_pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                                "chr_2", "pos_2", "R2", "Dprime")
    top_pairwise$dist <- abs(top_pairwise$pos_1 - top_pairwise$pos_2)
    top_pairwise_df[[i]] <- top_pairwise %>% mutate(R2 = round(R2, 3))

    top_diff_chr[[i]] <- filter(top_pairwise_df[[i]], chr_1 != chr_2) %>%
                         dplyr::select(R2) %>% unlist

  }; cat(" done.\n")

  all_top <- do.call("rbind", top_pairwise_df)
  #get sample size per pathway
  top_diff_num <- sapply(top_diff_chr, length)
  #get mean r2 value per pathway
  top_diff_mean <- sapply(top_diff_chr, mean)

  cat(sprintf("Calculated LD for %i total SNP pairs.\n", nrow(all_top)))
  cat(sprintf("%i SNP-SNP interactions in high-conf pathway %i (mean = %g)\n",
              top_diff_num, seq(top_bed), top_diff_mean))
  cat(sprintf("%i total interchromosomal SNP-SNP pairs.\n", sum(top_diff_num)))

  #remove original data objects to clear memory
  rm(top_list, top_ld_matrix)

  #============================================================================#
  # Permute random samples from original PLINK genotype data and calculate LD
  message("\n-------RANDOMLY SELECTED SNPS-------\n")

  # Large vector, time intensive
  start.time <- Sys.time()
  test <- read.plink(plinkF)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  rep.num <- 50L      #how many permutations to run
  sample.num <- 40L   #number of SNPs to pick for each permutation

  # Generating LD r2 stats for 500 permutations of 100 SNPs each
  # later will plot mean null r2 distribution via ggplot
  #null <- replicate(rep.num, {
  #  shuffle <- ld(test$genotypes[, sample(ncol(test$genotypes),
  #                              sample.num, replace=F)],
  #                stats="R.squared",
  #                depth=sample.num-1)
  #})

  null <- list()
  null_r2_matrix <- list()
  null_dprime_matrix <- list()
  null_pairwise_df <- list()
  null_diff_chr <- list()

  Sys.sleep(5)
  for (i in 1:rep.num) {
    cat(sprintf("Calculating LD within random sample matrix %i...", i))

    null[[i]] <- ld(test$genotypes[, sample(ncol(test$genotypes),
                                  sample.num, replace=F)],
                    stats=c("D.prime", "R.squared"),
                    depth=sample.num-1)

    null_r2_matrix[[i]] <- null[[i]]$R.squared
    null_dprime_matrix[[i]] <- null[[i]]$D.prime

    # Create dataframe containing pairwise distance calculations for each
    # LD SNP pair
    snp_map <- test$map

    # Turn each LD matrix into a data frame
    null_r2 <- as.matrix(null_r2_matrix[[i]]) #convert sparseMatrix to regular matrix
    null_r2 <- subset(melt(null_r2), value!=0) #melt df and remove '0's
    colnames(null_r2)[3] <- "R2"

    null_dprime <- as.matrix(null_dprime_matrix[[i]])
    null_dprime <- subset(melt(null_dprime), value!=0)
    colnames(null_dprime)[3] <- "Dprime"

    # Combine R2 and Dprime stats for each SNP-SNP pair
    null_stats <- merge(null_r2, null_dprime, by=c("Var1", "Var2"))

    # Generate pariwise distance table for each SNP-SNP pair
    colnames(null_stats)[1] <- "snp.name"
    snp_map <- subset(snp_map, select=c("snp.name", "chromosome", "position"))
    null_pairwise <- merge(snp_map, null_stats, by="snp.name")
    colnames(null_pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
    null_pairwise <- merge(snp_map, null_pairwise, by="snp.name")
    colnames(null_pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                                 "chr_2", "pos_2", "R2", "Dprime")
    # Calculate distance between SNP pairs
    null_pairwise$dist <- abs(null_pairwise$pos_1 - null_pairwise$pos_2)
    # Round r2 value to 3 decimal points
    null_pairwise_df[[i]] <- null_pairwise %>% mutate(R2 = round(R2, 3))
    null_pairwise_df[[i]] <- null_pairwise %>% mutate(Dprime = round(Dprime, 3))

    # Used to build null distruibution of mean R2 for SNP-SNP pairs on diff
    # chromosomes per sample 'pathway'
    null_diff_chr[[i]] <- filter(null_pairwise_df[[i]], chr_1 != chr_2) %>%
                          dplyr::select(R2) %>% unlist
    cat(" done.\n")
  }

  all_null <- do.call("rbind", null_pairwise_df)
  cat(sprintf("Calculated %i randomly permuted SNP interactions.\n",
              nrow(all_null)))

  null_diff_num <- sapply(null_diff_chr, length)
  #get mean r2 value per permutation
  null_diff_mean <- sapply(null_diff_chr, mean)
  rm(test)

  #============================================================================#
  ## PLOT STATS
  message("\n-------PLOTS-------\n")
  cat("Linkage disequilbrium statistic plots...")

  ####### high confidence
  #high confidence pairwise SNPs
  all_top$pathway_group <- "highconf"
  sub_top <- subset(all_top, R2 > 0.2)
  #filter based on matching/non-matching chromosome pairs
  all_top_same_chr <- filter(all_top, chr_1 == chr_2)
  all_top_diff_chr <- filter(all_top, chr_1 != chr_2)
  sub_top_same_chr <- filter(sub_top, chr_1 == chr_2)
  #write out tables
  write.table(all_top_same_chr,
              file=sprintf("%s/all_top_same_chr.txt", outDir),
                           col=T, row=F, quote=F, sep="\t")
  write.table(all_top_diff_chr,
              file=sprintf("%s/all_top_diff_chr.txt", outDir),
                           col=T, row=F, quote=F, sep="\t")

  ####### random
  #randomly selected pairwise SNPs
  #mean R2 per sample random 'pathway' to build null distribution
  #each row corresponds to each 'pathway' (nrow(null_diff_mean)==rep.num)
  all_null$pathway_group <- "null"
  sub_null <- subset(all_null, R2 > 0.2)
  #filter based on matching/non-matching chromosome pairs
  all_null_same_chr <- filter(all_null, chr_1 == chr_2)
  all_null_diff_chr <- filter(all_null, chr_1 != chr_2)
  sub_null_same_chr <- filter(sub_null, chr_1 == chr_2)

  ####### all
  #combine high-conf and null dataframes containing pairwise R2 values
  all_same_chr <- rbind(all_top_same_chr, all_null_same_chr)
  all_diff_chr <- rbind(all_top_diff_chr, all_null_diff_chr)
  sub_same_chr <- rbind(sub_top_same_chr, sub_null_same_chr)
  #divide distance axes to make them human-readable
  all_same_chr$dist <- all_same_chr$dist/1000000 ##divide by 1000000 (bp -> mb)
  all_diff_chr$dist <- all_diff_chr$dist/1000000
  sub_same_chr$dist <- sub_same_chr$dist/1000000

  ##############################################################################
  title_same <- paste("Degree of corrrelation per intrachromosomal SNP-SNP",
                      "\npair within the high-confidence pathway group vs.",
                      "\nrandomly selected SNP groups")
  title_diff <- paste("Degree of co-selection per interchromosomal SNP-SNP",
                      "\ninteraction within the high-confidence pathway",
                      "\ngroup vs. randomly selected SNP groups")

  ## for use with stat_summary(fun.data=box.style); allows white median line to
  ## appear after colouring and filling boxplots
  box.style <- function(x){
      return(c(y=median(x), ymin=median(x), ymax=median(x)))
    }

  ## for use with stat.summary(fun.data=give.n); displays sample size (N)
  ## courtesy of Bangyou at Stack Overflow
  give.n <- function(x){
      return(c(y = median(x)*1.30, label = length(x)))
      # experiment with the multiplier to find the perfect position
    }

  # boxplot of r2 per high confidence pathway
  blah <- melt(top_diff_chr)
  ggplot(blah, aes(x=factor(L1), y=value, color="#F8766D", fill="#F8766D")) +
       geom_boxplot() +
       stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
                    fun.data=box.style) +
       stat_summary(geom="text", color="white", fun.data=give.n,
                    position=position_dodge(width=0.75)) +
       scale_y_continuous(name=bquote("Pairwise LD value (" *r^2*")")) +
       scale_x_discrete(name="# of interchromosomal SNP-SNP pairs per pathway") +
       ggtitle(paste("Degree of co-selection per interchromosomal SNP-SNP",
                     "\npair within the high-confidence pathways")) +
       theme(axis.text.x=element_text(vjust=0.4, hjust=1)) +
       theme_set(theme_minimal()) +
       theme(plot.title=element_text(hjust=0.5),
             text=element_text(size=17),
             legend.position="none",
             panel.grid.major.x=element_blank()) +
       geom_hline(yintercept=0.1, colour="grey", linetype="dashed", size=1)
  ggsave("high_interchr_bar.png", width=11)

  #################### PLOT 1: Null distribution vs real #######################
  # Plotting mean R2 value for random SNP-SNP pairs on diff chromosome
  # against real mean R2 value for high-conf SNPs
  # interchromosomal SNP pairs used as proxy for co-selection/genetic ixn

  null <- null_diff_mean # mean r2 values per permutation (n=500)
  real <- mean(all_top_diff_chr$R2) # mean r2 for all high-conf pathways (n=1)

  ggplot(as.data.frame(null), aes(x=as.data.frame(null))) +
      geom_histogram(aes(y=..density..),
                         colour="white", fill="grey", bins=15) +
      geom_density(alpha=0.3, colour="#00BFC4", fill="#00BFC4") +
      scale_x_continuous(name=bquote("Pairwise LD value (mean " *r^2*")")) +
      scale_y_continuous(name="Density") +
      ggtitle(title_diff) +
      #xintercept = mean R2 for high-conf SNP-SNP pairs on diff chrs (real R2)
      geom_vline(aes(xintercept=real),
                     color="#F8766D", linetype="dashed", size=1) +
      theme_set(theme_minimal()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank())
  ggsave("dist_null_real.png", width=8, height=7)

  #calculate pvalue by quantifying all permuted mean r2 values greater than
  #the real mean r2 values, divided by the total number of replicates
  pval <- function(null, real) {
    perm.p = ( length(which(null > mean(real)))+1 ) / (length(null)+1)
    return(perm.p)
  }

  perm.p = pval(null, real)
  cat(sprintf("P-value for permuted sample vs. real mean r2 = %g\n", perm.p))

  ########################## PLOT 2: R2 boxplots ###############################
  # R2 boxplot without distance bins
  dat <- sub_same_chr
  ggplot(dat, aes(x=pathway_group, y=-log10(R2))) +
        geom_boxplot(outlier.colour=NULL, aes(colour=pathway_group,
                     fill=pathway_group)) +
        stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
                     fun.data=box.style) +
        scale_y_continuous(name=bquote("Pairwise LD value (-log10("*r^2*"))")) +
        #scale_y_continuous(name=bquote("Pairwise LD value ("*r^2*")")) +
        ggtitle(title_same) +
        theme_set(theme_minimal()) +
        theme(plot.title=element_text(hjust=0.5),
              text=element_text(size=17),
              legend.position="top",
              legend.title=element_blank(),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_blank()) +
        scale_x_discrete(labels=paste("N=", table(dat$pathway_group), sep=""))
    ggsave("sub_same_log.png", width=7.5, height=7.5)

  ######################### PLOT 3: Density plots ##############################
  # Get mean R2 for each pathway group and plot density
  ##NOTE:only do this for intrachromosomal SNPs (all_same_chr + sub_same_chr)
  meanvals <- ddply(dat, "pathway_group", summarise, r2_mean=mean(R2))
  ggplot(dat, aes(x=R2)) +
      geom_density(aes(group=pathway_group, colour=pathway_group,
                       fill=pathway_group), alpha=0.3) +
      scale_x_continuous(name=bquote("Pairwise LD value (" *r^2*")")) +
      scale_y_continuous(name="Density") +
      ggtitle(title_same) +
      theme_set(theme_minimal()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank()) +
      geom_vline(data=meanvals, aes(xintercept=r2_mean,
                 colour=pathway_group), linetype="dashed", size=1)
  ggsave("sub_density.png", width=8, height=7)

  # Get significance values
  dat <- all_diff_chr
  wilcox <- wilcox.test(R2 ~ pathway_group, data=dat, alternative="greater")
  ttest <- t.test(dat$R2 ~ pathway_group, var.equal=T)

  print(wilcox)
  print(ttest)

  ######################## PLOT 5: Linear regression ###########################
  # Plot linear regression of R2 vs distance for both pathway groups
  ##NOTE:only regressing intrachromosomal SNPs (all_same_chr + sub_same_chr)
  ggplot(dat, aes(x=dist, y=R2, colour=pathway_group)) +
        geom_point() +
        geom_smooth(method=lm, se=T) +
        scale_y_continuous(name=bquote("Pairwise LD value (" *r^2*")")) +
        scale_x_continuous(name="Chromosomal distance (Mb)") +
        ggtitle(title_same) +
        theme_set(theme_minimal()) +
        theme(plot.title=element_text(hjust=0.5),
              text=element_text(size=17),
              legend.position="top",
              legend.title=element_blank(),
              panel.grid.major.x=element_blank()) +
       coord_cartesian(xlim=c(0, 0.25), ylim=c(0,1)) #zoom in
   ggsave("sub_regr_zoom.png", width=8, height=7)

  # Determine significance of regression
  fit <- lm(formula=R2 ~ dist + pathway_group, data=dat)
  summary(fit)

#==============================================================================
  # Plots no longer used
  # R2 boxplot by distance (first cut distance into 4 equally sized groups)
#  (q<-quantile(dat$dist, seq(0, 1, 0.25)))
#  dat$dist_ranges <- cut(dat$dist, breaks=q, include.lowest=T)

#  ggplot(dat, aes(x=dist_ranges, y=-log10(R2), colour=pathway_group,
#               fill=pathway_group)) +
#      geom_boxplot(outlier.colour=NULL) +
#      stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
#                   fun.data=box.style) +
#      scale_x_discrete("Chromosomal distance (Mb)") +
#      scale_y_continuous(name=bquote("-log10("*R^2*")")) +
#      ggtitle(title_same) +
#      theme_set(theme_minimal()) +
#      theme(plot.title=element_text(hjust=0.5), text=element_text(size=17),
#            legend.position="top", legend.title=element_blank(),
#            panel.grid.major.x=element_blank())

  ##contour plots - clearer visaul of bimodality
#  ggplot(sub_same_chr, aes(dist, R2)) +
#        geom_density_2d() +
#        geom_point() +
#        scale_y_continuous(name=expression(paste("R"^"2"))) +
#        scale_x_continuous(name="Distance (bp)") +
#        ggtitle(paste("Distribution of pairwise R2 measure as a function of",
#                      "distance\n within the high-confidence pathway group")) +
#         geom_hline(aes(yintercept=0.2), colour="red") +
#        theme(plot.title=element_text(hjust=0.5), text=element_text(size=17))

  #Plot R2 for all high-conf pathways
  ##fill vectors with NA before melting list
#  max_length <- max(unlist(lapply(top_ld_r2, length)))
#  r2_filled <- lapply(top_ld_r2, function(x){ans <- rep(NA, length=max_length);
#                                             ans[1:length(x)]<- x;
#                                             return(ans)})

#  r2_filled <- do.call(cbind, r2_filled)
#  r2_melt <- melt(r2_filled)

  ##rename columns for better plot visualization
#  colnames(r2_melt) <- c("SNP_pairs", "pathway", "R2")
#  r2_melt <- na.omit(r2_melt)
#  r2_melt$pathway <- paste("High-conf_path", r2_melt$pathway, sep=" ")

#  ggplot(r2_melt, aes(x=SNP_pairs, y=R2)) +
#       facet_wrap(~pathway, nrow=3, ncol=6) +
#       scale_y_continuous(name=expression(paste("R"^"2"))) +
#       scale_x_continuous(name="Number of SNP-SNP pairs in pathway") +
#       geom_point(aes(colour=pathway)) +
#       geom_hline(yintercept=0.5, colour="red") +
#       theme_bw() +
#       theme(text=element_text(size=15))
  #ggsave()
}
