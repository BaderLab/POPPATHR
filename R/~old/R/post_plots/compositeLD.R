require(SNPRelate)
require(reshape2)
require(data.table)
require(dplyr)
require(tools)
require(ggplot2)
require(RColorBrewer)
require(cowplot)

pop1 <- "CEU"
pop2 <- "ASW"
popNames <- c("European", "African")
hcInDir <- "hc_snps"
lcInDir <- "lc_snps"

source("../../../bin/R/themePublication.R")

statDir <- sprintf("res_%s_%s_compLD_absD", hcInDir, lcInDir)
if (!file.exists(statDir)) dir.create(statDir)

LDstats <- function(hcInDir, lcInDir, popNames=NULL, statDir) {

  calcLD <- function(snpDir, popcode) {

    bed2gds <- c(); genos <- c(); ld <- c()
    pairwise.df <- list(); diff.df <- list()

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

    # Convert PLINK files to GDS objects
    # NOTE: order is important!
    for (i in 1:length(bim)) {
      cat(sprintf("\n*Reading genotype information for the %s pathway\n",
        basename(file_path_sans_ext(bim[i]))))
      Sys.sleep(0)
      bed2gds[[i]] <- snpgdsBED2GDS(bed[i], fam[i], bim[i],
                                    sprintf("%s/pathway_%i", statDir, i))

      # Open GDS object per pathway
      genos[[i]] <- snpgdsOpen(bed2gds[i])

      # Subset genotypes by population
      if (popcode == 0) {

        # Measure LD (use sample.id to specify samples to calculate LD with)
        cat("\n*Measuring LD per each pairwise SNP-SNP combination\n")
        ld[[i]] <- snpgdsLDMat(genos[[i]], slide=-1, method="composite")

      } else {

        famF <- fread(fam[[i]], h=F, data.table=F)
        pop <- filter(famF, V6 == popcode) %>% select(V2) %>% unlist

        cat(sprintf("\n*Keeping %i genotypes from %s population(s)\n",
            length(pop),
            if (popcode == 1) {pop1}
            else if (popcode == 2) {pop2} ))

        # Measure LD (use sample.id to specify samples to calculate LD with)
        cat("*Measuring LD per each pairwise SNP-SNP combination\n")
        ld[[i]] <- snpgdsLDMat(genos[[i]], sample.id=pop,
                              slide=-1, method="corr")

      }

      # reshape matrix and get relevant values (ie remove duplicates and
      # LD calculated for the same snp-snp interaction)
      ld.mat <- ld[[i]]$LD
      snp.id <- ld[[i]]$snp.id
      colnames(ld.mat) <- snp.id
      rownames(ld.mat) <- snp.id
      ld.mat[lower.tri(ld.mat)] <- NA  # get the upper triangular matrix
      ld.mat[!upper.tri(ld.mat)] <- NA # remove the diagonal

      # Melt matrix to a data frame
      pairwise <- melt(ld.mat)
      pairwise <- na.omit(pairwise)
      colnames(pairwise) <- c("snp_1", "snp_2", "corr")

      # SNP positions (use fread to specify columns to read in)
      snp.map <- fread(bim[[i]], h=F, select=c(1:2),
                       col.names=c("chr_1", "snp_1"),
                       data.table=F)

      # Assign chromosomal position to each SNP
      pairwise <- merge(pairwise, snp.map, by="snp_1")
      colnames(snp.map) <- c("chr_2", "snp_2")
      pairwise <- merge(pairwise, snp.map, by="snp_2")
      pairwise$pathway <- sprintf("pathway_%i", i)

      # Round the correlation value
      pairwise.df[[i]] <- pairwise %>% mutate(corr=round(corr, 3))

      # Keep only inter-chr interactions
      cat("\n*Subsetting for inter-chromosomal SNP-SNP pairs only\n")
      diff.df[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2)
      cat(sprintf("\tTotal # of inter-chromosomal SNP-SNP interactions: %i\n",
          nrow(diff.df[[i]])))

      # Close gds files and unlink
      snpgdsClose(genos[[i]])
      unlink(sprintf("%s/pathway_%i", statDir, i))

      cat(" done.\n\n")
    }

    all.pairs  <<- do.call("rbind", pairwise.df)
    diff.pairs <<- do.call("rbind", diff.df)
    diff.pairs$corr <<- abs(diff.pairs$corr)
    diff.num   <<- sapply(diff.df, nrow) # all inter-chr SNP-SNP pairs per path

    cat(sprintf("Finished inter-chr LD analysis for %i pathways.\n", length(bim)))
    cat(sprintf("Calculated LD for %i total SNP pairs.\n", nrow(all.pairs)))
    cat(sprintf(" --> %i total interchromosomal pairs.\n\n", sum(diff.num)))

  }

  # Set alternate population names if provided
  if (!is.null(popNames)) {
    pop1.name <- popNames[1]
    pop2.name <- popNames[2]
  } else {
    pop1.name <- pop1
    pop2.name <- pop2
  }

  sink(sprintf("%s/interChrLDanalysis.log", statDir))

  # Calculate inter-chromosomal LD stats for confidently enriched pathways
  cat("=======================================================================")
  cat(paste("\n*Calculating inter-chromosomal LD between the confidently",
            "enriched pathway variants...\n"))

  # all population genotypes
  cat("========== ALL POPULATIONS ==========\n")
  calcLD(snpDir=hcInDir, popcode=0)
  hc.diff.pairs <- diff.pairs
  hc.diff.pairs$set <- "Enriched"
  hc.diff.num <- diff.num
  cat(" done.\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop1))
  Sys.sleep(3)
  calcLD(snpDir=hcInDir, popcode=1)
  hc.diff.pairs.pop1 <- diff.pairs
  hc.diff.pairs.pop1$set <- "Enriched"
  hc.diff.pairs.pop1$pop <- pop1.name
  hc.diff.num.pop1 <- diff.num
  cat(" done.\n")

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop2))
  Sys.sleep(3)
  calcLD(snpDir=hcInDir, popcode=2)
  hc.diff.pairs.pop2 <- diff.pairs
  hc.diff.pairs.pop2$set <- "Enriched"
  hc.diff.pairs.pop2$pop <- pop2.name
  hc.diff.num.pop2 <- diff.num
  cat(" done.\n")

  # Calculate inter-chromosomal LD stats for nonenriched pathways
  cat("=======================================================================")
  cat(paste("\n*Calculating inter-chromosomal LD between the nonenriched",
            "pathway variants...\n"))

  # all population genotypes
  cat("\n========== ALL POPULATIONS ==========\n")
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=0)
  lc.diff.pairs <- diff.pairs
  lc.diff.pairs$set <- "Nonenriched"
  lc.diff.num <- diff.num
  cat(" done.\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop1))
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=1)
  lc.diff.pairs.pop1 <- diff.pairs
  lc.diff.pairs.pop1$set <- "Nonenriched"
  lc.diff.pairs.pop1$pop <- pop1.name
  lc.diff.num.pop1 <- diff.num
  cat(" done.\n")

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ONLY ==========\n", pop2))
  Sys.sleep(3)
  calcLD(snpDir=lcInDir, popcode=2)
  lc.diff.pairs.pop2 <- diff.pairs
  lc.diff.pairs.pop2$set <- "Nonenriched"
  lc.diff.pairs.pop2$pop <- pop2.name
  lc.diff.num.pop2 <- diff.num
  cat(" done.\n")
  cat("=======================================================================")

  # Write out tables of interactions per pathway
  path_names <- list.files(path=hcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  hc.ixns <- data.frame(pathway=path_names,
                        hc_ixns=hc.diff.num,
                        hc_ixns_pop1=hc.diff.num.pop1,
                        hc_ixns_pop2=hc.diff.num.pop2)
  write.table(hc.ixns, sprintf("%s/hc_num_interactions_pathway.txt", statDir),
              col.names=T, row.names=F, quote=F, sep="\t")

  path_names <- list.files(path=lcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  lc.ixns <- data.frame(pathway=path_names,
                        lc_ixns=lc.diff.num,
                        lc_ixns_pop1=lc.diff.num.pop1,
                        lc_ixns_pop2=lc.diff.num.pop2)
  write.table(lc.ixns, sprintf("%s/lc_num_interactions_pathway.txt", statDir),
              col.names=T, row.names=F, quote=F, sep="\t")

  ## PLOT STATS
  cat("\n*Generating inter-chromosomal LD analysis plots.\n")
  ############################### PLOT 1 #####################################
  # Pairwise inter r2 of all enriched against all nonenriched pathways
  # Set common variables
  title <- "Enrichment of SNP-SNP interactions within the ancestry-enriched pathways"
  r2.xaxis.title <- bquote(bold("LD value per SNP-SNP pair"))

  dat <- rbind(hc.diff.pairs, lc.diff.pairs)

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=corr, colour=set, fill=set)) +
            geom_density(alpha=0.2) +
            xlab(r2.xaxis.title) +
            ylab(bquote(bold("Density"))) +
            scale_colour_Publication() +
            scale_fill_Publication()

  # 1b) eCDF plot (cumulative density at each r2)
  p2 <- ggplot(dat, aes(x=corr, colour=set)) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_Publication() +
          scale_fill_Publication()

  cat("*Generating plot 1...")
  both_1  <- plot_grid(p1, p2, labels=c("A", "B"), ncol=2)
  title_1 <- ggdraw() + draw_label(title, fontface="bold")
  both_1  <- plot_grid(title_1, both_1, ncol=1, rel_heights=c(0.1, 1))
  filename_1 <- sprintf("%s/density_ecdf_hc_lc_r2.png", statDir)
  save_plot(filename_1, both_1, base_height=3.5, base_width=10.5,
            base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename_1))

  # Calculate significance via KS test
  pval <- ks.test(hc.diff.pairs$corr,
                  lc.diff.pairs$corr,
                  alternative="less")
  cat(sprintf("\n\t**p-value via KS test (less)= %g**\n", pval$p.value))

  ############################### PLOT 2 #####################################
  # Pairwise inter r2 per enriched and nonenriched pathway separately
  title <- paste("Pathway-specific enrichment of",
                 "inter-chromosomal SNP-SNP interactions")

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$pathway)))

  # Function to determine significance of inter-chrom LD per each
  # enriched pathway vs. cumulative set of nonenriched pathways
  # via the KS test (alternative=less{the CDF of x lies above that of y})
  getPvals <- function(hcSet, lcSet) {
    for (i in 1:length(unique(hcSet$pathway))) {
      path_ld[[i]] <<- subset(hcSet, pathway==sprintf("pathway_%s", i))
      ks_pvals[[i]] <<- ks.test(path_ld[[i]]$corr,
                                lcSet$corr,
                                alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Write out each p value alongside enriched pathway names
  path_names <- list.files(path=hcInDir, pattern="*.snps$", full.names=F)
  path_names <- gsub("\\%.*", "", path_names)
  path_names <- gsub("_", " ", path_names)

  # 2a) Density distribution plot
  p3 <- ggplot(dat, aes(x=corr, colour=pathway, fill=pathway)) +
          facet_grid(set ~ .) +
          geom_density(alpha=0.2) +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Density"))) +
          scale_fill_manual(values=npal) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text = element_text(face="bold"))

  # 2b) eCDF plot (cumulative density at each r2)
  p4 <- ggplot(dat, aes(x=corr, colour=pathway)) +
          facet_grid(set ~ .) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text = element_text(face="bold"))

  # 2c,d) Density and eCDF at x-axis > 0.2
  p5 <- p3 + xlim(0.5, 1)
  p6 <- p4 + xlim(0.5, 1)

  cat("\n*Generating plot 2...")
  both_2 <- plot_grid(p3, p4, p5, p6, labels=c("A", "B", "C", "D"),
                      ncol=2, nrow=2)
  title_2 <- ggdraw() + draw_label(title, fontface="bold")
  both_2 <- plot_grid(title_2, both_2, ncol=1, rel_heights=c(0.1, 1))
  filename_2 <- sprintf("%s/inter_den_ecdf_hc_lc.png", statDir)
  save_plot(filename_2, both_2, base_height=8, base_width=9.5,
            base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename_2))

  # Calculate p values
  path_ld <- list()
  ks_pvals <- list()

  getPvals(hcSet=hc.diff.pairs, lcSet=lc.diff.pairs)
  pvals <- unlist(ks_pvals)
  path_pvals <- data.frame(pathway=path_names,
                           pval=as.data.frame(pvals),
                           bonf=p.adjust(pvals, method="bonferroni"),
                           fdr=p.adjust(pvals, method="BH"))

   path_pvals <- path_pvals[order(path_pvals$fdr),]
   filename_p1 <- sprintf("%s/hc_pvals_per_pathway_alt-l.txt", statDir)
   write.table(format(path_pvals, digits=3), file=filename_p1,
               col=T, row=F, quote=F, sep="\t")
   cat(sprintf("*Table of pathway p-values written to %s.\n", filename_p1))

  ############################### PLOT 3 #####################################
  # Pairwise inter r2 per enriched and nonenriched pathway separately
  # and stratified by population
  title <- paste("Pathway-specific enrichment of",
              "inter-chromosomal SNP-SNP interactions per population")

  dat <- rbind(hc.diff.pairs.pop1, lc.diff.pairs.pop1,
              hc.diff.pairs.pop2, lc.diff.pairs.pop2)

  # 3a) Density distribution plot
  p7 <- ggplot(dat, aes(x=corr, colour=pathway, fill=pathway)) +
          facet_grid(set ~ pop) +
          geom_density(alpha=0.2) +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Density"))) +
          scale_fill_manual(values=npal) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

  # 3b) eCDF plot (cumulative density at each r2)
  p8 <- ggplot(dat, aes(x=corr, colour=pathway)) +
           facet_grid(set ~ pop) +
           stat_ecdf() +
           xlab(r2.xaxis.title) +
           ylab(bquote(bold("Cumulative density"))) +
           scale_colour_manual(values=npal) +
           theme(legend.position="none",
                 strip.text=element_text(face="bold"))

   # 3c,d) Density and eCDF at x-axis > 0.2
   p9  <- p7 + xlim(0.5, 1)
   p10 <- p8 + xlim(0.5, 1)

   cat("\n*Generating plot 3...")
   both_3 <- plot_grid(p7, p8, p9, p10, labels=c("A", "B", "C", "D"),
                       ncol=2, nrow=2)
   title_3 <- ggdraw() + draw_label(title, fontface='bold')
   both_3  <- plot_grid(title_3, both_3, ncol=1, rel_heights=c(0.1, 1))
   filename_3 <- sprintf("%s/inter_pop-strat_den_ecdf_hc_lc.png", statDir)
   save_plot(filename_3, both_3, base_height=10, base_width=13.5,
             base_aspect_ratio=1.2)
   cat(sprintf(" saved to %s.\n", filename_2))

   # Calculate p values per population
   # Pop 1 calculations
   path_ld <- list()
   ks_pvals <- list()

   getPvals(hc=hc.diff.pairs.pop1, lc=lc.diff.pairs.pop1)
   pvals <- unlist(ks_pvals)
   path_pvals_pop1 <- data.frame(pathway=path_names,
                                 pval=as.data.frame(pvals),
                                 bonf=p.adjust(pvals, method="bonferroni"),
                                 fdr=p.adjust(pvals, method="BH"))

   path_pvals_pop1 <- path_pvals_pop1[order(path_pvals_pop1$fdr),]
   filename_p2 <- sprintf("%s/hc_pvals_per_pathway_alt-l_%s.txt", statDir, pop1)
   write.table(format(path_pvals_pop1, digits=3), file=filename_p2,
               col=T, row=F, quote=F, sep="\t")
   cat(sprintf("*Table of pathway p-values written to %s.\n", filename_p2))

   # Pop 2 calculations
   path_ld <- list()
   ks_pvals <- list()

   getPvals(hc=hc.diff.pairs.pop2, lc=lc.diff.pairs.pop2)
   pvals <- unlist(ks_pvals)
   path_pvals_pop2 <- data.frame(pathway=path_names,
                                 pval=as.data.frame(pvals),
                                 bonf=p.adjust(pvals, method="bonferroni"),
                                 fdr=p.adjust(pvals, method="BH"))

   path_pvals_pop2 <- path_pvals_pop2[order(path_pvals_pop2$fdr),]
   filename_p3 <- sprintf("%s/hc_pvals_per_pathway_alt-l_%s.txt", statDir, pop2)
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
   filename_p4 <- sprintf("%s/hc_pvals_merge_%s_%s.txt", statDir, pop1, pop2)
   write.table(format(pval_merge, digits=3), file=filename_p4,
               col=T, row=F, quote=F, sep="\t")
   cat(sprintf("Merged p-value tables written to %s.\n", filename_p4))

   sink()
}
