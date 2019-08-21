#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#'
#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param lcInDir (char) path to files with low-confidence pathway SNP lists
#' @param snp2geneF (char) path to file with snp2gene mappings. Output of
#' 		mapSNP2gene() (found in GWAS2Pathway)
#' @param outDir (char) path to output directory
#'
#' @return none
#' @export
#'
LDstatsBPM <- function(hcInDir, lcInDir, snp2geneF, outDir) {

# Read in snp2gene.txt file to assign SNPs to genes
snp2gene <- read.table(snp2geneF, h=F, as.is=T)
snp2gene <- snp2gene[,-3]

  calcLD <- function(snpDir, popcode) {

    ss1 <- c(); ss2 <- c()
    snp.map1 <- c(); snp.map2 <- c()
    ld.calc <- list(); pairwise.df <- list(); diff <- list()
    diff.r2 <- list(); diff.dp <- list()

    # Read PLINK files for each pathway
    bed <- list.files(path=snpDir, pattern="*.bed", full.names=T)
    bim <- list.files(path=snpDir, pattern="*.bim", full.names=T)
    fam <- list.files(path=snpDir, pattern="*.fam", full.names=T)

    # Matrices of all possible pathway x pathway combinations
    bed.pair <- t(combn(bed, 2))
    bim.pair <- t(combn(bim, 2))
    fam.pair <- t(combn(fam, 2))

    for (i in 1:nrow(bed.pair)) {
      cat(sprintf("#### BETWEEN-PATHWAY INTERACTION %i ####\n", i))
      cat(sprintf("*Reading pathway %s x pathway %s PLINK sets\n",
            basename(file_path_sans_ext(bed.pair[i,1])),
            basename(file_path_sans_ext(bed.pair[i,2]))))

      ss1[[i]] <- read.plink(bed.pair[i,1], bim.pair[i,1], fam.pair[i,1])
      ss2[[i]] <- read.plink(bed.pair[i,2], bim.pair[i,2], fam.pair[i,2])

      # Subset genotypes by population
      if (popcode == 0) {
        # For first pathway of pathwayxpathway interaction
        pop <- which(ss1[[i]]$fam$affected != popcode)
        ss1[[i]]$genotypes <- ss1[[i]]$genotypes[pop, ]
        cat(sprintf("\n*Keeping %i individuals in pathway one of interaction\n",
                length(pop)))
        print(ss1[[i]]$genotypes)

        # For second pathway of pathwayxpathway interaction
        pop <- which(ss2[[i]]$fam$affected != popcode)
        ss2[[i]]$genotypes <- ss2[[i]]$genotypes[pop, ]
        cat(sprintf("\n*Keeping %i individuals in pathway two of interaction\n",
                length(pop)))
        print(ss2[[i]]$genotypes)

      } else {
        # For first pathway of pathwayxpathway interaction
        pop <- which(ss1[[i]]$fam$affected == popcode)
        ss1[[i]]$genotypes <- ss1[[i]]$genotypes[pop, ]
        cat(sprintf("\n*Keeping %i %s genotypes in pathway one of interaction\n",
            length(pop),
            if (popcode == 1) {pop1}
            else if (popcode == 2) {pop2} ))
        print(ss1[[i]]$genotypes)

        # For second pathway of pathwayxpathway interaction
        pop <- which(ss2[[i]]$fam$affected == popcode)
        ss2[[i]]$genotypes <- ss2[[i]]$genotypes[pop, ]
        cat(sprintf("\n*Keeping %i %s genotypes in pathway two of interaction\n",
            length(pop),
            if (popcode == 1) {pop1}
            else if (popcode == 2) {pop2} ))
        print(ss2[[i]]$genotypes)
      }

      cat("*Calculating pathway x pathway LD statistics...\n")
      ld.calc[[i]] <- ld(x=ss1[[i]]$genotypes, y=ss2[[i]]$genotypes,
                         stats=c("D.prime", "R.squared"))

      snp.map1 <- ss1[[i]]$map #genomic location of each SNP
      snp.map2 <- ss2[[i]]$map

      r2 <- as.matrix(ld.calc[[i]]$R.squared) #convert sparseMatrix to regular matrix
      r2 <- melt(r2) #melt matrix to data frame
      colnames(r2)[3] <- "R.squared"

      # Create dataframe containing pairwise distance calculations for each
      # SNP-SNP pair
      dp <- as.matrix(ld.calc[[i]]$D.prime)
      dp <- melt(dp)
      colnames(dp)[3] <- "D.prime"

      all.stats <- merge(r2, dp, by=c("Var1", "Var2"))
      colnames(all.stats)[1:2] <- c("snp.name.1", "snp.name.2")

      snp.map1 <- subset(snp.map1, select=c("snp.name", "chromosome", "position"))
      colnames(snp.map1)[1] <- "snp.name.1"
      snp.map2 <- subset(snp.map2, select=c("snp.name", "chromosome", "position"))
      colnames(snp.map2)[1] <- "snp.name.2"

      pairwise <- merge(snp.map1, all.stats, by="snp.name.1")
      colnames(pairwise)[1:3] <- c("snp_1", "chr_1", "pos_1")

      pairwise <- merge(snp.map2, pairwise, by="snp.name.2")
      colnames(pairwise) <- c("snp_2", "chr_2", "pos_2", "snp_1",
                              "chr_1", "pos_1", "R.squared", "D.prime")
      pairwise <- pairwise[,c(4:6, 1:3, 7, 8)]

      pairwise$pathway_pair1 <- basename(file_path_sans_ext(bed.pair[i,1]))
      pairwise$pathway_pair2 <- basename(file_path_sans_ext(bed.pair[i,2]))
      pairwise$ixn_num <- sprintf("interaction_%i", i)

      # Assign gene to each SNP
      colnames(snp2gene) <- c("snp_1", "gene_1")
      pairwise.df <- merge(pairwise, snp2gene, by="snp_1")
      colnames(snp2gene) <- c("snp_2", "gene_2")
      pairwise.df <- merge(pairwise.df, snp2gene, by="snp_2")
      cat(sprintf("\tTotal SNP-SNP interactions: %i\n", nrow(pairwise.df)))

      cat("*Removing any matching genes in pathway x pathway interaction...\n")
      remove <- intersect(pairwise.df$gene_1, pairwise.df$gene_2)
      no.match <- pairwise.df[!(pairwise.df$gene_1 %in% remove),]
      no.match2 <- no.match[!(no.match$gene_2 %in% remove),]

      cat("*Filtering for inter-chromosomal interactions...\n")
      diff[[i]] <- filter(no.match2, chr_1 != chr_2)
      cat(sprintf("\tNew number of interactions: %i\n", nrow(diff[[i]])))

      diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
      diff.dp[[i]] <- select(diff[[i]], D.prime) %>% unlist
      cat("done.\n\n")
    }

    diff.pairs <<- do.call("rbind", diff)
    diff.num <<- sapply(diff.r2, length)   #all SNP-SNP pairs per interaction

    cat(sprintf("Finished inter-chr LD analysis for %i pathway x pathway interactions.\n",
        nrow(bed.pair)))

  }

  if (!is.null(popNames)) {
    pop1.name <- popNames[1]
    pop2.name <- popNames[2]
  } else {
    pop1.name <- pop1
    pop2.name <- pop2
  }

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
  # Get dataframe listing each unique pathway x pathway interaction
  pathway_ixns <- unique(hc.diff.pairs[,c("pathway_pair1", "pathway_pair2")])

  hc.ixns <- data.frame(interaction=pathway_ixns,
                        hc_ixns=hc.diff.num,
                        hc_ixns_pop1=hc.diff.num.pop1,
                        hc_ixns_pop2=hc.diff.num.pop2)
  write.table(hc.ixns, sprintf("%s/hc_num_interactions_pathway-pathway.txt",
              outDir), col.names=T, row.names=F, quote=F, sep="\t")

  pathway_ixns <- unique(lc.diff.pairs[,c("pathway_pair1", "pathway_pair2")])

  lc.ixns <- data.frame(interaction=pathway_ixns,
                        lc_ixns=lc.diff.num,
                        lc_ixns_pop1=lc.diff.num.pop1,
                        lc_ixns_pop2=lc.diff.num.pop2)
  write.table(lc.ixns, sprintf("%s/lc_num_interactions_pathway-pathway.txt",
              outDir), col.names=T, row.names=F, quote=F, sep="\t")

  lc.diff.num[lc.diff.num==0] <- NA
  lc.diff.num <- na.omit(lc.diff.num)

  lc.diff.num.pop1[lc.diff.num.pop1==0] <- NA
  lc.diff.num.pop1 <- na.omit(lc.diff.num.pop1)

  lc.diff.num.pop2[lc.diff.num.pop2==0] <- NA
  lc.diff.num.pop2 <- na.omit(lc.diff.num.pop2)

  ## PLOT STATS
  cat("\n*Generating inter-chromosomal LD analysis plots.\n")
  ############################### PLOT 1 #####################################
  # Pairwise inter r2 of all enriched against all unenriched pathways
  # Set common variables
  title <- "Enrichment of SNP-SNP interactions between the selection-enriched pathways"
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
  title <- paste("Between pathway enrichment of",
                 "inter-chromosomal SNP-SNP interactions")

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$ixn_num)))

  # Function to determine significance of inter-chrom LD per selection-enriched
  # between-pathway interaction vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  getPvals <- function(dat, pop.name) {

    # Separate df into 'Enriched' and 'Unenriched' interactions
    if (is.null(pop.name) == TRUE) {
      enriched <- filter(dat, set=="Enriched")
      unenriched <- filter(dat, set=="Unenriched")
    } else { #population-stratified
      enriched <- filter(dat, pop==pop.name & set=="Enriched")
      unenriched <- filter(dat, pop==pop.name & set=="Unenriched")
   }

    for (i in 1:length(unique(enriched$ixn_num))) {
      # Subset for each enriched pathway
      enrich_ixn_ld[[i]] <<- subset(enriched, ixn_num==sprintf("interaction_%s", i))
      # Calculate KS pvalue for each enriched pathway against the entire
      # set of unenriched pathways
      ks_pvals[[i]] <<- ks.test(enrich_ixn_ld[[i]]$TEST.STAT,
                                unenriched$TEST.STAT,
                                alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Get dataframe listing each unique pathway x pathway interaction
  pathway_ixns <- unique(hc.diff.pairs[,c("pathway_pair1", "pathway_pair2")])

  # 2a) Density distribution plot
  p3 <- ggplot(dat, aes(x=TEST.STAT, colour=ixn_num, fill=ixn_num)) +
          facet_grid(set ~ .) +
          geom_density(alpha=0.2) +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Density"))) +
          scale_fill_manual(values=npal) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

  # 2b) eCDF plot (cumulative density at each r2)
  p4 <- ggplot(dat, aes(x=TEST.STAT, colour=ixn_num)) +
          facet_grid(set ~ .) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

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
  enrich_ixn_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop=NULL)
  pvals <- unlist(ks_pvals)
  ixn_pvals <- data.frame(interaction=pathway_ixns,
                          pval=as.data.frame(pvals),
                          bonf=p.adjust(pvals, method="bonferroni"),
                          fdr=p.adjust(pvals, method="BH"))

  ixn_pvals <- ixn_pvals[order(ixn_pvals$fdr),]
  filename_p1 <- sprintf("%s/hc_pvals_per_interaction_alt-l.txt", outDir)
  write.table(format(ixn_pvals, digits=3), file=filename_p1,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of interaction p-values written to %s.\n", filename_p1))

  ############################### PLOT 3 #####################################
  # Pairwise inter r2 per enriched and unenriched pathway separately
  # and stratified by population
  title <- paste("Between pathway enrichment of",
                 "inter-chromosomal SNP-SNP interactions per population")

  dat <- rbind(hc.diff.pairs.pop1, lc.diff.pairs.pop1,
               hc.diff.pairs.pop2, lc.diff.pairs.pop2)

  # Rename dataframe column based on chosen LD statistic, R.squared or D.prime
  # in order to ensure consistency calling the correct column
  names(dat) <- gsub(statistic, "TEST.STAT", names(dat))
  cat(sprintf("*Generating results based on '%s' statistic\n", statistic))

  # 3a) Density distribution plot
  p7 <- ggplot(dat, aes(x=TEST.STAT, colour=ixn_num, fill=ixn_num)) +
         facet_grid(set ~ pop) +
         geom_density(alpha=0.2) +
         xlab(r2.xaxis.title) +
         ylab(bquote(bold("Density"))) +
         scale_fill_manual(values=npal) +
         scale_colour_manual(values=npal) +
         theme(legend.position="none",
               strip.text=element_text(face="bold"))

  # 3b) eCDF plot (cumulative density at each r2)
  p8 <- ggplot(dat, aes(x=TEST.STAT, colour=ixn_num)) +
          facet_grid(set ~ pop) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

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
  enrich_ixn_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop.name=pop1.name)
  pvals <- unlist(ks_pvals)
  ixn_pvals_pop1 <- data.frame(interaction=pathway_ixns,
                               pval=as.data.frame(pvals),
                               bonf=p.adjust(pvals, method="bonferroni"),
                               fdr=p.adjust(pvals, method="BH"))

  ixn_pvals_pop1 <- ixn_pvals_pop1[order(ixn_pvals_pop1$fdr),]
  filename_p2 <- sprintf("%s/hc_pvals_per_interaction_alt-l_%s.txt", outDir, pop1)
  write.table(format(ixn_pvals_pop1, digits=3), file=filename_p2,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of interaction p-values written to %s.\n", filename_p2))

  # Pop 2 calculations
  enrich_ixn_ld <- list()
  ks_pvals <- list()

  getPvals(dat, pop.name=pop2.name)
  pvals <- unlist(ks_pvals)
  ixn_pvals_pop2 <- data.frame(interaction=pathway_ixns,
                               pval=as.data.frame(pvals),
                               bonf=p.adjust(pvals, method="bonferroni"),
                               fdr=p.adjust(pvals, method="BH"))

  ixn_pvals_pop2 <- ixn_pvals_pop2[order(ixn_pvals_pop2$fdr),]
  filename_p3 <- sprintf("%s/hc_pvals_per_interaction_alt-l_%s.txt", outDir, pop2)
  write.table(format(ixn_pvals_pop2, digits=3), file=filename_p3,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("*Table of interaction p-values written to %s.\n", filename_p3))

  # Merge both p value dataframes
  colnames(ixn_pvals_pop1)[3:5] <- paste(colnames(ixn_pvals_pop1)[3:5],
                                          sep="_", pop1)
  colnames(ixn_pvals_pop2)[3:5] <- paste(colnames(ixn_pvals_pop2)[3:5],
                                          sep="_", pop2)
  pval_merge <- merge(ixn_pvals_pop1, ixn_pvals_pop2,
                      by=c("interaction.pathway_pair1",
                           "interaction.pathway_pair2"))

  pval_merge <- pval_merge[order(pval_merge[,5]), ]
  filename_p4 <- sprintf("%s/hc_pvals_merge_%s_%s.txt", outDir, pop1, pop2)
  write.table(format(pval_merge, digits=3), file=filename_p4,
              col=T, row=F, quote=F, sep="\t")
  cat(sprintf("Merged p-value tables written to %s.\n", filename_p4))
}
