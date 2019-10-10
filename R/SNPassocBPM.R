#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#'
#' @param enrichDir (char) path to selection-enriched pathway SNP lists
#' @param unenrichDir (char) path to unenriched pathway SNP lists
#' @param pop1 (char) character code for the first population (controls).
#' @param pop2 (char) character code for the second population (cases).
#' @param snp2geneF (char) path to file with snp2gene mappings. Output of
#' 		mapSNP2gene() (found in GWAS2Pathway)
#' @param outDir (char) path to output directory
#'
#' @return none
#' @export
#'

SNPassocBPM <- function(enrichDir, unenrichDir, pop1, pop2,
                        snp2geneF, outDir) {
  # Read in snp2gene file to assign SNPs to genes to remove any
  # matching genes in a given pathway-pathway pair
  snp2gene <- read.table(snp2geneF, h=FALSE, as.is=TRUE)
  snp2gene <- snp2gene[,-3]

  calcLD <- function(snpDir, popcode) {

    ss1 <- c(); ss2 <- c()
    snp.map1 <- c(); snp.map2 <- c()
    ld.calc <- list(); pairwise.df <- list(); diff <- list()
    diff.r2 <- list(); diff.dp <- list()

    # Read PLINK files for each pathway
    bed <- list.files(path=snpDir, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snpDir, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snpDir, pattern="*.fam", full.names=TRUE)

    # Matrices of all possible pathway x pathway combinations
    bed.pair <- t(combn(bed, 2))
    bim.pair <- t(combn(bim, 2))
    fam.pair <- t(combn(fam, 2))

    for (i in 1:nrow(bed.pair)) {
      path.name1 <- substr(basename(bim.pair[i,1]), 0, nchar(basename(bim.pair[i,1]))-4)
      path.name2 <- substr(basename(bim.pair[i,2]), 0, nchar(basename(bim.pair[i,2]))-4)

      cat(sprintf("BETWEEN-PATHWAY INTERACTION %i (%s total)\n", i, nrow(bim.pair)))
      cat(sprintf("*Reading SNPs in pathway %s x pathway %s\n", path.name1, path.name2))

      ss1[[i]] <- read.plink(bed.pair[i,1], bim.pair[i,1], fam.pair[i,1])
      ss2[[i]] <- read.plink(bed.pair[i,2], bim.pair[i,2], fam.pair[i,2])

      # Subset genotypes by population
      # For first pathway of pathway x pathway interaction
      pop <- which(ss1[[i]]$fam$affected == popcode)
      ss1[[i]]$genotypes <- ss1[[i]]$genotypes[pop, ]
      cat(sprintf("\n%i %s genotypes in pathway %s\n",
          length(pop), ifelse(popcode == 1, pop1, pop2), path.name1))
      print(ss1[[i]]$genotypes)

      # For second pathway of pathway x pathway interaction
      pop <- which(ss2[[i]]$fam$affected == popcode)
      ss2[[i]]$genotypes <- ss2[[i]]$genotypes[pop, ]
      cat(sprintf("\n%i %s genotypes in pathway %s\n",
          length(pop), ifelse(popcode == 1, pop1, pop2), path.name2))
      print(ss2[[i]]$genotypes)

      cat("\n*Calculating pathway x pathway LD statistics...\n")
      ld.calc[[i]] <- ld(x=ss1[[i]]$genotypes, y=ss2[[i]]$genotypes,
                         stats="R.squared")

      snp.map1 <- ss1[[i]]$map # genomic location of each SNP
      snp.map2 <- ss2[[i]]$map

      r2 <- as.matrix(ld.calc[[i]]) # convert sparseMatrix to regular matrix
      r2 <- melt(r2) # melt matrix to data frame
      colnames(r2)[3] <- "R.squared"

      # Create dataframe with pairwise distance calculations for each SNP-SNP pair
      colnames(r2)[1:2] <- c("snp.name.1", "snp.name.2")

      snp.map1 <- subset(snp.map1, select=c("snp.name", "chromosome", "position"))
      colnames(snp.map1)[1] <- "snp.name.1"
      snp.map2 <- subset(snp.map2, select=c("snp.name", "chromosome", "position"))
      colnames(snp.map2)[1] <- "snp.name.2"

      pairwise <- join(snp.map1, r2, by="snp.name.1")
      colnames(pairwise)[1:3] <- c("snp_1", "chr_1", "pos_1")
      pairwise <- join(snp.map2, pairwise, by="snp.name.2")
      colnames(pairwise) <- c("snp_2", "chr_2", "pos_2", "snp_1",
                              "chr_1", "pos_1", "R.squared")
      pairwise <- pairwise[,c(4:6, 1:3, 7)]

      pairwise$pathway_pair1 <- path.name1
      pairwise$pathway_pair2 <- path.name2
      pairwise$ixn_num <- sprintf("interaction_%i", i)

      # Assign SNPs to genes
      colnames(snp2gene) <- c("snp_1", "gene_1")
      pairwise.df <- join(pairwise, snp2gene, by="snp_1")
      colnames(snp2gene) <- c("snp_2", "gene_2")
      pairwise.df <- join(pairwise.df, snp2gene, by="snp_2")
      cat(sprintf("\tTotal SNP-SNP interactions: %i\n", nrow(pairwise.df)))

      cat("*Removing any matching genes in pathway x pathway interaction...\n")
      remove <- intersect(pairwise.df$gene_1, pairwise.df$gene_2)
      no.match <- pairwise.df[!(pairwise.df$gene_1 %in% remove),]
      no.match2 <- no.match[!(no.match$gene_2 %in% remove),]

      cat("*Filtering for inter-chromosomal interactions...\n")
      diff[[i]] <- filter(no.match2, chr_1 != chr_2)
      cat(sprintf("\tNew number of interactions: %i\n", nrow(diff[[i]])))

      diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
      cat("done.\n\n")
    }

    diff.pairs <<- do.call("rbind", diff)
    diff.num <<- sapply(diff.r2, length)  # all SNP-SNP pairs per interaction

    cat(sprintf("Finished inter-chr LD analysis for %i pathway x pathway interactions.\n",
        nrow(bed.pair)))
  }

  # Calculate inter-chromosomal LD stats between enriched pathways
  cat("=======================================================================")
  cat("**Measuring trans-chromosomal LD between selection-enriched pathways\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop1))
  calcLD(snpDir=enrichDir, popcode=1)
  enrich.pop1 <- diff.pairs
  enrich.pop1$set <- "Enriched"
  enrich.pop1$pop <- pop1
  enrich.num.pop1 <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop2))
  calcLD(snpDir=enrichDir, popcode=2)
  enrich.pop2 <- diff.pairs
  enrich.pop2$set <- "Enriched"
  enrich.pop2$pop <- pop2
  enrich.num.pop2 <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # Calculate inter-chromosomal LD stats between unenriched pathways
  cat("=======================================================================")
  cat("**Measuring trans-chromosomal LD between unenriched pathways\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop1))
  calcLD(snpDir=unenrichDir, popcode=1)
  unenrich.pop1 <- diff.pairs
  unenrich.pop1$set <- "Unenriched"
  unenrich.pop11$pop <- pop1
  unenrich.num.pop1 <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop2))
  calcLD(snpDir=unenrichDir, popcode=2)
  unenrich.pop2 <- diff.pairs
  unenrich.pop2$set <- "Unenriched"
  unenrich.pop2$pop <- pop2
  unenrich.num.pop2 <- diff.num
  cat(" done.\n")
  Sys.sleep(3)
  cat("=======================================================================")

  # Write out tables of SNP-SNP interactions per pathway-pathway combo
  col_pop1 <- sprintf("num_ixns_%s", pop1)
  col_pop2 <- sprintf("num_ixns_%s", pop2)

  ## Enriched
  ixn_enrich <- unique(enrich.pop1[,c("pathway_pair1", "pathway_pair2")])
  enrich_df <- data.frame(interaction=ixn_enrich,
                          pop1=enrich.num.pop1,
                          pop2=enrich.num.pop2)
  colnames(enrich_df)[2] <- col_pop1
  colnames(enrich_df)[3] <- col_pop2

  write.table(enrich_df,
    sprintf("%s/enrich_num_interactions_pathway-pathway.txt", outDir),
      col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  ixn_unenrich <- unique(unenrich.pop1[,c("pathway_pair1", "pathway_pair2")])
  unenrich_df <- data.frame(interaction=ixn_unenrich,
                            pop1=unenrich.num.pop1,
                            pop2=unenrich.num.pop2)
  colnames(unenrich_df)[2] <- col_pop1
  colnames(unenrich_df)[3] <- col_pop2

  write.table(unenrich_df,
    sprintf("%s/unenrich_num_interactions_pathway-pathway.txt", outDir),
      col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  unenrich.num.pop1[unenrich.num.pop1==0] <- NA
  unenrich.num.pop1 <- na.omit(unenrich.num.pop1)

  unenrich.num.pop2[unenrich.num.pop2==0] <- NA
  unenrich.num.pop2 <- na.omit(unenrich.num.pop2)

  ## PLOT STATS
  cat("\n*Generating inter-chromosomal LD plots and analyses.\n")
  # Function to determine significance of inter-chrom LD per selection-enriched
  # between-pathway interaction vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  getPvals <- function(dat, pop_name) {
    # Separate df into 'Enriched' and 'Unenriched' interactions
    if (is.null(pop.name) == TRUE) {
      enriched <- filter(dat, set=="Enriched")
      unenriched <- filter(dat, set=="Unenriched")
    } else { #population-stratified
      enriched <- filter(dat, pop==pop_name & set=="Enriched")
      unenriched <- filter(dat, pop==pop_name & set=="Unenriched")
   }

    for (i in 1:length(unique(enriched$ixn_num))) {
      # Subset for each enriched pathway
      enrich_ixn_ld[[i]] <<- subset(enriched, ixn_num==sprintf("interaction_%s", i))
      # Calculate KS pvalue for each enriched pathway against the entire
      # set of unenriched pathways
      ks_pvals[[i]] <<- ks.test(enrich_ixn_ld[[i]]$R.squared,
                                unenriched$R.squared,
                                alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Population-specific pairwise inter-chromosomal r2 per enriched pathway
  title <- "Pathway-specific SNP-SNP coevolution between enriched pathways"
  r2.xaxis.title <- "LD value per SNP-SNP pair"

  # Merge all results together
  dat <- rbind(enrich.pop1, unenrich.pop1, enrich.pop2, unenrich.pop2)

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$ixn_num)))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=R.squared, colour=ixn_num, fill=ixn_num)) +
         facet_grid(set ~ pop) +
         geom_density(alpha=0.2) +
         xlab(r2.xaxis.title) +
         ylab(bquote(bold("Density"))) +
         scale_fill_manual(values=npal) +
         scale_colour_manual(values=npal) +
         theme(legend.position="none",
               strip.text=element_text(face="bold"))

  # 1b) eCDF plot (cumulative density at each r2)
  p2 <- ggplot(dat, aes(x=R.squared, colour=ixn_num)) +
          facet_grid(set ~ pop) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

  # 3c,d) Density and eCDF at x-axis > 0.2
  p3 <- p1 + xlim(0.2, 1)
  p4 <- p2 + xlim(0.2, 1)

  cat("\n*Generating plot...")
  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
  title <- ggdraw() + draw_label(title, fontface="bold")
  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
  filename <- sprintf("%s/between_snp-snp_interchrom_ld_r2.png", outDir)
  save_plot(filename, both, base_height=10, base_width=13.5, base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename))

  # Calculate p values per population
  enrich_path_ld <- list()
  ks_pvals <- list()
  res <- list()

  cat("*Determining significant SNP-SNP coevolution within enriched pathways.\n")
  calcCoev <- function(pop) {
    # p value function defined above
    getPvals(dat=dat, pop_name=pop)

    # Generate p value dataframe
    pvals <- unlist(ks_pvals)
    pvals_df <- data.frame(interaction=ixn_enrich,
                           pval=as.data.frame(pvals),
                           fdr=p.adjust(pvals, method="BH"))

    pvals_df <- pvals_df[order(pvals_df$fdr),]
    sig_paths <- filter(pvals_df, fdr <= 0.2)
    cat(sprintf("Pathway-pathway interactions with significant coevolution in %s (FDR<=0.2, N=%s):\n%s\n",
        pop, nrow(sig_paths), paste(sig_paths$pathway, collapse="\n")))
    colnames(pvals_df)[2:3] <- paste(colnames(pvals_df)[2:3], sep="_", pop)
    res[[pop]] <<- pvals_df

    # Write out results
    filename_1 <- sprintf("%s/enrich_coevolution_pval_%s.txt", outDir, pop)
    write.table(format(pvals_df, digits=3), file=filename_1,
        col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    cat(sprintf("*Table of pathway p-values written to %s.\n", filename_1))
  }
  # Run for both populations
  mapply(calcCoev, c(pop1, pop2))

  # Merge both p value dataframes
  pval_merge <- join(res[[1]], res[[2]], by=c("interaction.pathway_pair1", "interaction.pathway_pair2"))
  filename_2 <- sprintf("%s/enrich_coevolution_pval_merge_%s_%s.txt", outDir, pop1, pop2)
  write.table(format(pval_merge, digits=3), file=filename_2,
      col=TRUE, row=FALSE, quote=FALSE, sep="\t")
  cat(sprintf("\n*Merged p-value tables written to %s.\n", filename_2))
}
