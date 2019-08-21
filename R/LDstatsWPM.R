#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#' https://bioconductor.org/packages/release/bioc/manuals/snpStats/man/snpStats.pdf

#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param lcInDir (char) path to files with low-confidence pathway SNP lists
#' @param statistic (char) function will generate both R squared and D prime
#'    LD statistics; choose either 'R.squared' or 'D.prime' for plotting and
#'    p value calculation (default="R.squared")
#'    NOTE: syntax is important!
#' @param popNames (char) optional - character vector to set alternate
#'    population names in plots eg. c("European", "African") rather than
#'    "CEU" and "ASW"
#' @param outDir (char) path to output directory
#'
#' @return none
#' @export
#'
LDstatsWPM <- function(hcInDir, lcInDir, statistic="R.squared",
                    popNames=NULL, outDir) {

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

  # Skip LD calculation for pathways with <= 1 SNP
  skip <- sapply(bim, readLines)
  skip <- which(sapply(skip, length) <= 1)

  if (length(skip) > 0) {
    cat(sprintf("**Skipping LD calculation for pathway %s-not enough SNPs\n\n",
        skip))
    bed <- bed[-skip]
    bim <- bim[-skip]
    fam <- fam[-skip]
  } else {
    cat("All pathways have enough SNPs for LD calculation. Commencing...\n\n")
  }

  # Convert PLINK files to snpStats input format
  # Output object is a list with 3 elements ($genotypes, $fam, $map)
  # NOTE: order is important!
  for (i in 1:length(bim)) {
    cat(sprintf("*Reading pathway %s PLINK set\n",
          basename(file_path_sans_ext(bed[i]))))

    ss[[i]] <- read.plink(bed[i], bim[i], fam[i])

    # Subset genotypes by population
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

    cat("*Calculating LD statistics\n")
    ld.calc[[i]] <- ld(ss[[i]]$genotypes,
                       stats=c("D.prime", "R.squared"),
                       depth=ncol(ss[[i]]$genotypes)-1)
    snp.map <- ss[[i]]$map #genomic location of each SNP

    # Turn each LD matrix into a data frame
    r2 <- as.matrix(ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
    lowerTriangle(r2, byrow=F, diag=T) <- NA #replace lower triangle and diag with NA
                                             #NOTE:updated 03/23/2018
    r2 <- na.omit(melt(r2)) #keep all non-NA values
    colnames(r2)[3] <- "R.squared"

    # Create dataframe containing pairwise distance calculations for each
    # SNP-SNP pair
    dp <- as.matrix(ld.calc[[i]]$D.prime)
    lowerTriangle(dp, byrow=F, diag=T) <- NA #replace lower triangle and diag with NA
                                             #NOTE:updated 03/23/2018
    dp <- na.omit(melt(dp)) #keep all non-NA values
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

    # NOTE: rounding may affect p value caluclation by KS test (03/08/2018)
    #pairwise.df[[i]] <- pairwise %>% mutate(R.squared=round(R.squared, 3))
    #pairwise.df[[i]] <- pairwise %>% mutate(D.prime=round(D.prime, 3))

    pairwise.df[[i]] <- pairwise %>% mutate(R.squared=R.squared)
    pairwise.df[[i]] <- pairwise %>% mutate(D.prime=D.prime)

    cat("*Subsetting for interchromosomal SNP-SNP pairs\n")
    diff[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2)

    diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
    diff.dp[[i]] <- select(diff[[i]], D.prime) %>% unlist
    cat("done.\n\n")
  }

  all.pairs <<- do.call("rbind", pairwise.df)
  diff.pairs <<- do.call("rbind", diff)
  diff.r2.mean <<- sapply(diff.r2, mean) #mean r2 value per pathway
  diff.num <<- sapply(diff.r2, length)   #all SNP-SNP pairs per path
  #diff.dp.mean <<- sapply(diff.dp, mean) #mean dprime per pathway

  cat(sprintf("Finished inter-chr LD analysis for %i pathways.\n", length(bim)))
  cat(sprintf("Calculated LD for %i total SNP pairs.\n", nrow(all.pairs)))
  cat(sprintf(" --> %i total interchromosomal pairs.\n\n", sum(diff.num)))
}

  # How all pairwise combinations are calculated
  # n is the number of SNPs in each pathway, divided by 2 to account for
  # duplicate pairs ("n choose 2")
  # combos <- function(n) {
  #   return( n*(n-1)/2 ) }

  # Set alternate population names if provided
  if (!is.null(popNames)) {
    pop1.name <- popNames[1]
    pop2.name <- popNames[2]
  } else {
    pop1.name <- pop1
    pop2.name <- pop2
  }

  sink(sprintf("%s/interChrLDanalysis.log", outDir))

  # Calculate inter-chromosomal LD stats for confidently enriched pathways
  cat("=======================================================================")
  cat(paste("\n*Calculating inter-chromosomal LD between the selection-",
            "enriched pathway variants...\n"))

  # all population genotypes
  cat("========== ALL POPULATIONS ==========\n")
  calcLD(snpDir=hcInDir, popcode=0)
  hc.diff.pairs <- diff.pairs
  hc.diff.pairs$set <- "Enriched"
  hc.diff.num <- diff.num
  cat(" done.\n")
  saveRDS(hc.diff.pairs, sprintf("%s/hc.diff.pairs.rds", outDir))

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop1))
  Sys.sleep(3)
  calcLD(snpDir=hcInDir, popcode=1)
  hc.diff.pairs.pop1 <- diff.pairs
  hc.diff.pairs.pop1$set <- "Enriched"
  hc.diff.pairs.pop1$pop <- pop1.name
  hc.diff.num.pop1 <- diff.num
  cat(" done.\n")
  saveRDS(hc.diff.pairs.pop1, sprintf("%s/hc.diff.pairs.pop1.rds", outDir))

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop2))
  Sys.sleep(3)
  calcLD(snpDir=hcInDir, popcode=2)
  hc.diff.pairs.pop2 <- diff.pairs
  hc.diff.pairs.pop2$set <- "Enriched"
  hc.diff.pairs.pop2$pop <- pop2.name
  hc.diff.num.pop2 <- diff.num
  cat(" done.\n")
  saveRDS(hc.diff.pairs.pop2, sprintf("%s/hc.diff.pairs.pop2.rds", outDir))

  # Calculate inter-chromosomal LD stats for unenriched pathways
  cat("=======================================================================")
  cat(paste("\n*Calculating inter-chromosomal LD between the unenriched",
            "pathway variants...\n"))

  # all population genotypes
  cat("\n========== ALL POPULATIONS ==========\n")
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=0)
  lc.diff.pairs <- diff.pairs
  lc.diff.pairs$set <- "Unenriched"
  lc.diff.num <- diff.num
  cat(" done.\n")
  saveRDS(lc.diff.pairs, sprintf("%s/lc.diff.pairs.rds", outDir))

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop1))
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=1)
  lc.diff.pairs.pop1 <- diff.pairs
  lc.diff.pairs.pop1$set <- "Unenriched"
  lc.diff.pairs.pop1$pop <- pop1.name
  lc.diff.num.pop1 <- diff.num
  cat(" done.\n")
  saveRDS(lc.diff.pairs.pop1, sprintf("%s/lc.diff.pairs.pop1.rds", outDir))

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop2))
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=2)
  lc.diff.pairs.pop2 <- diff.pairs
  lc.diff.pairs.pop2$set <- "Unenriched"
  lc.diff.pairs.pop2$pop <- pop2.name
  lc.diff.num.pop2 <- diff.num
  cat(" done.\n")
  saveRDS(lc.diff.pairs.pop2, sprintf("%s/lc.diff.pairs.pop2.rds", outDir))
  cat("=======================================================================")

  # Write out tables of interactions per pathway
  path_names <- list.files(path=hcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  hc.ixns <- data.frame(pathway=path_names,
                        hc_ixns=hc.diff.num,
                        hc_ixns_pop1=hc.diff.num.pop1,
                        hc_ixns_pop2=hc.diff.num.pop2)
  write.table(hc.ixns, sprintf("%s/hc_num_interactions_pathway.txt", outDir),
              col.names=T, row.names=F, quote=F, sep="\t")

  path_names <- list.files(path=lcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  lc.ixns <- data.frame(pathway=path_names,
                        lc_ixns=lc.diff.num,
                        lc_ixns_pop1=lc.diff.num.pop1,
                        lc_ixns_pop2=lc.diff.num.pop2)
  write.table(lc.ixns, sprintf("%s/lc_num_interactions_pathway.txt", outDir),
              col.names=T, row.names=F, quote=F, sep="\t")

  ## PLOT STATS
  cat("\n*Generating inter-chromosomal LD analysis plots.\n")
  ############################### PLOT 1 #####################################
  # Pairwise inter r2 of all enriched against all unenriched pathways
  # Set common variables
  title <- "Enrichment of SNP-SNP interactions within the selection-enriched pathways"
  r2.xaxis.title <- bquote(bold("LD value per SNP-SNP pair"))

  dat <- rbind(hc.diff.pairs, lc.diff.pairs)

  # Rename dataframe column based on chosen LD statistic, R.squared or D.prime
  # in order to ensure consistency calling the correct column
  names(dat) <- gsub(statistic, "TEST.STAT", names(dat))
  cat(sprintf("*Generating results based on '%s' statistic\n", statistic))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=TEST.STAT, colour=set, fill=set)) +
            geom_density(alpha=0.2) +
            xlab(r2.xaxis.title) +
            ylab(bquote(bold("Density"))) +
            scale_colour_Publication(guide=FALSE) +
            scale_fill_Publication(guide=FALSE)

  # 1b) eCDF plot (cumulative density at each r2)
  p2 <- ggplot(dat, aes(x=TEST.STAT, colour=set)) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_Publication(guide=FALSE) +
          scale_fill_Publication(guide=FALSE)

  cat("*Generating plot 1...")
  both_1  <- plot_grid(p1, p2, labels=c("A", "B"), ncol=2)
  title_1 <- ggdraw() + draw_label(title, fontface="bold")
  both_1  <- plot_grid(title_1, both_1, ncol=1, rel_heights=c(0.1, 1))
  filename_1 <- sprintf("%s/density_ecdf_hc_lc_%s.tiff", outDir, statistic)
  save_plot(filename_1, both_1, base_height=3.5, base_width=10.5,
            base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename_1))

  # Calculate significance via KS test
  pval <- ks.test(filter(dat, set=="Enriched") %>% select(TEST.STAT) %>% unlist,
                  filter(dat, set=="Unenriched") %>% select(TEST.STAT) %>% unlist,
                  alternative="less")
  cat(sprintf("\n\t**p-value via KS test (less)= %g**\n", pval$p.value))

  ############################### PLOT 2 #####################################
  # Pairwise inter r2 per enriched and unenriched pathway separately
  title <- paste("Pathway-specific enrichment of",
                 "inter-chromosomal SNP-SNP interactions")

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$pathway)))

  # Function to determine significance of inter-chrom LD per each
  # enriched pathway vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  getPvals <- function(dat, pop.name) {

    # Separate df into 'Enriched' and 'Unenriched' pathways
    if (is.null(pop.name) == TRUE) {
      enriched <- filter(dat, set=="Enriched")
      unenriched <- filter(dat, set=="Unenriched")
    } else { #population-stratified
      enriched <- filter(dat, pop==pop.name & set=="Enriched")
      unenriched <- filter(dat, pop==pop.name & set=="Unenriched")
   }

    for (i in 1:length(unique(enriched$pathway))) {
      # Subset for each enriched pathway
      enrich_path_ld[[i]] <<- subset(enriched, pathway==sprintf("pathway_%s", i))
      # Calculate KS pvalue for each enriched pathway against the entire
      # set of unenriched pathways
      ks_pvals[[i]] <<- ks.test(enrich_path_ld[[i]]$TEST.STAT,
                                unenriched$TEST.STAT,
                                alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Write out each p value alongside enriched pathway names
  path_names <- list.files(path=hcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  # 2a) Density distribution plot
  p3 <- ggplot(dat, aes(x=TEST.STAT, colour=pathway, fill=pathway)) +
          facet_grid(set ~ .) +
          geom_density(alpha=0.2) +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Density"))) +
          scale_fill_manual(values=npal) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text = element_text(face="bold"))

  # 2b) eCDF plot (cumulative density at each r2)
  p4 <- ggplot(dat, aes(x=TEST.STAT, colour=pathway)) +
          facet_grid(set ~ .) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text = element_text(face="bold"))

  # 2c,d) Density and eCDF at x-axis > 0.2
  p5 <- p3 + xlim(0.2, 1)
  p6 <- p4 + xlim(0.2, 1)

  cat("\n*Generating plot 2...")
  both_2 <- plot_grid(p3, p4, p5, p6, labels=c("A", "B", "C", "D"),
                      ncol=2, nrow=2)
  title_2 <- ggdraw() + draw_label(title, fontface="bold")
  both_2 <- plot_grid(title_2, both_2, ncol=1, rel_heights=c(0.1, 1))
  filename_2 <- sprintf("%s/inter_den_ecdf_hc_lc_%s.tiff", outDir, statistic)
  save_plot(filename_2, both_2, base_height=8, base_width=9.5,
            base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename_2))

  # Calculate p values
  enrich_path_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop=NULL)
  pvals <- unlist(ks_pvals)
  path_pvals <- data.frame(pathway=path_names,
                           pval=as.data.frame(pvals),
                           bonf=p.adjust(pvals, method="bonferroni"),
                           fdr=p.adjust(pvals, method="BH"))

  path_pvals <- path_pvals[order(path_pvals$fdr),]
  filename_p1 <- sprintf("%s/hc_pvals_per_pathway_alt-l.txt", outDir)
  write.table(format(path_pvals, digits=3), file=filename_p1,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of pathway p-values written to %s.\n", filename_p1))

  ############################### PLOT 3 #####################################
  # Pairwise inter r2 per enriched and unenriched pathway separately
  # and stratified by population
  title <- paste("Pathway-specific enrichment of",
                 "inter-chromosomal SNP-SNP interactions per population")

  dat <- rbind(hc.diff.pairs.pop1, lc.diff.pairs.pop1,
               hc.diff.pairs.pop2, lc.diff.pairs.pop2)

  # Rename dataframe column based on chosen LD statistic, R.squared or D.prime
  # in order to ensure consistency calling the correct column
  names(dat) <- gsub(statistic, "TEST.STAT", names(dat))
  cat(sprintf("*Generating results based on '%s' statistic\n", statistic))

  # 3a) Density distribution plot
  p7 <- ggplot(dat, aes(x=TEST.STAT, colour=pathway, fill=pathway)) +
         facet_grid(set ~ pop) +
         geom_density(alpha=0.2) +
         xlab(r2.xaxis.title) +
         ylab(bquote(bold("Density"))) +
         scale_fill_manual(values=npal) +
         scale_colour_manual(values=npal) +
         theme(legend.position="none",
               strip.text = element_text(face="bold"))

  # 3b) eCDF plot (cumulative density at each r2)
  p8 <- ggplot(dat, aes(x=TEST.STAT, colour=pathway)) +
          facet_grid(set ~ pop) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text = element_text(face="bold"))

  # 3c,d) Density and eCDF at x-axis > 0.2
  p9  <- p7 + xlim(0.2, 1)
  p10 <- p8 + xlim(0.2, 1)

  cat("\n*Generating plot 3...")
  both_3 <- plot_grid(p7, p8, p9, p10, labels=c("A", "B", "C", "D"),
                      ncol=2, nrow=2)
  title_3 <- ggdraw() + draw_label(title, fontface='bold')
  both_3  <- plot_grid(title_3, both_3, ncol=1, rel_heights=c(0.1, 1))
  filename_3 <- sprintf("%s/inter_pop-strat_den_ecdf_hc_lc_%s.tiff",
                    outDir, statistic)
  save_plot(filename_3, both_3, base_height=10, base_width=13.5,
            base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename_2))

  # Calculate p values per population
  # Pop 1 calculations
  enrich_path_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop.name=pop1.name)
  pvals <- unlist(ks_pvals)
  path_pvals_pop1 <- data.frame(pathway=path_names,
                                pval=as.data.frame(pvals),
                                bonf=p.adjust(pvals, method="bonferroni"),
                                fdr=p.adjust(pvals, method="BH"))

  path_pvals_pop1 <- path_pvals_pop1[order(path_pvals_pop1$fdr),]
  filename_p2 <- sprintf("%s/hc_pvals_per_pathway_alt-l_%s.txt", outDir, pop1)
  write.table(format(path_pvals_pop1, digits=3), file=filename_p2,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of pathway p-values written to %s.\n", filename_p2))

  # Pop 2 calculations
  enrich_path_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop.name=pop2.name)
  pvals <- unlist(ks_pvals)
  path_pvals_pop2 <- data.frame(pathway=path_names,
                                pval=as.data.frame(pvals),
                                bonf=p.adjust(pvals, method="bonferroni"),
                                fdr=p.adjust(pvals, method="BH"))

  path_pvals_pop2 <- path_pvals_pop2[order(path_pvals_pop2$fdr),]
  filename_p3 <- sprintf("%s/hc_pvals_per_pathway_alt-l_%s.txt", outDir, pop2)
  write.table(format(path_pvals_pop2, digits=3), file=filename_p3,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of pathway p-values written to %s.\n", filename_p3))

  # Merge both p value dataframes
  colnames(path_pvals_pop1)[2:4] <- paste(colnames(path_pvals_pop1)[2:4],
                                          sep="_", pop1)
  colnames(path_pvals_pop2)[2:4] <- paste(colnames(path_pvals_pop2)[2:4],
                                          sep="_", pop2)
  pval_merge <- merge(path_pvals_pop1, path_pvals_pop2, by="pathway")

  pval_merge <- pval_merge[order(pval_merge[,2]), ]
  filename_p4 <- sprintf("%s/hc_pvals_merge_%s_%s.txt", outDir, pop1, pop2)
  write.table(format(pval_merge, digits=3), file=filename_p4,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("Merged p-value tables written to %s.\n", filename_p4))

  sink()
}
