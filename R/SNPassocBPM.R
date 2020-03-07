#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#'
#' @param pop_one (char) character code for the first population.
#' @param pop_two (char) character code for the second population.
#' @param snp2gene_file (char) path to file with snp2gene mappings
#'    (output of SNP2gene.R).
#' @param ASSOC_FDR (interger) FDR value to define significant pathway coevolution.
#' @param enrich_folder (char) path to selection-enriched pathway SNP lists
#' @param unenrich_folder (char) path to unenriched pathway SNP lists
#' @param output_folder (char) path to output directory
#'
#' @return none
#' @export
#'

SNPassocBPM <- function(pop_one, pop_two, snp2gene_file, ASSOC_FDR,
                        enrich_folder, unenrich_folder,
                        output_folder) {
  # Read in snp2gene file to assign SNPs to genes to remove any
  # matching genes in a given pathway-pathway pair
  snp2gene <- read.table(snp2gene_file, h=FALSE, as.is=TRUE)
  snp2gene <- snp2gene[,-3]

  calcLD <- function(snp_folder, popcode) {

    ss1 <- c(); ss2 <- c()
    snp.map1 <- c(); snp.map2 <- c()
    ld.calc <- list(); pairwise.df <- list(); diff <- list()
    diff.r2 <- list(); diff.dp <- list()

    # Read PLINK files for each pathway
    bed <- list.files(path=snp_folder, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snp_folder, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snp_folder, pattern="*.fam", full.names=TRUE)

    # Matrices of all possible pathway x pathway combinations
    bed.pair <- t(combn(bed, 2))
    bim.pair <- t(combn(bim, 2))
    fam.pair <- t(combn(fam, 2))

    for (i in 1:nrow(bed.pair)) {
      path.name1 <- substr(basename(bim.pair[i,1]), 0, nchar(basename(bim.pair[i,1]))-4)
      path.name2 <- substr(basename(bim.pair[i,2]), 0, nchar(basename(bim.pair[i,2]))-4)

      cat(sprintf("BETWEEN-PATHWAY INTERACTION %i (%s total)\n", i, nrow(bim.pair)))
      cat(sprintf("* Reading SNPs in pathway %s x pathway %s\n", path.name1, path.name2))

      ss1[[i]] <- read.plink(bed.pair[i,1], bim.pair[i,1], fam.pair[i,1])
      ss2[[i]] <- read.plink(bed.pair[i,2], bim.pair[i,2], fam.pair[i,2])

      # Subset genotypes by population
      # For first pathway of pathway x pathway interaction
      pop <- which(ss1[[i]]$fam$affected == popcode)
      ss1[[i]]$genotypes <- ss1[[i]]$genotypes[pop, ]
      cat(sprintf("\n%i %s genotypes in pathway %s\n", length(pop), ifelse(popcode == 1, pop_one, pop_two), path.name1))
      print(ss1[[i]]$genotypes)

      # For second pathway of pathway x pathway interaction
      pop <- which(ss2[[i]]$fam$affected == popcode)
      ss2[[i]]$genotypes <- ss2[[i]]$genotypes[pop, ]
      cat(sprintf("\n%i %s genotypes in pathway %s\n",
          length(pop), ifelse(popcode == 1, pop_one, pop_two), path.name2))
      print(ss2[[i]]$genotypes)

      cat("\n* Calculating pathway x pathway SNP association statistics...\n")
      ld.calc[[i]] <- ld(x=ss1[[i]]$genotypes, y=ss2[[i]]$genotypes, stats="R.squared")

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
      colnames(pairwise) <- c("snp_2", "chr_2", "pos_2", "snp_1", "chr_1", "pos_1", "R.squared")
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

      cat("*Filtering for trans-chromosomal interactions...\n")
      diff[[i]] <- filter(no.match2, chr_1 != chr_2)
      cat(sprintf("\tNew number of interactions: %i\n", nrow(diff[[i]])))

      diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
      cat("done.\n\n")
    }

    diff.pairs <<- do.call("rbind", diff)
    diff.num <<- sapply(diff.r2, length)  # all SNP-SNP pairs per pathway interaction

    cat(sprintf("Finished inter-chr LD analysis for %i pathway x pathway interactions.\n",
        nrow(bed.pair)))
  }

  # Calculate trans-chromosomal LD stats between enriched pathways
  cat("* Calculating trans-chromosomal SNP association between selection-enriched pathways\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop_one))
  calcLD(snp_folder=enrich_folder, popcode=1)
  enrich.pop_one <- diff.pairs
  enrich.pop_one$set <- "Enriched"
  enrich.pop_one$pop <- pop_one
  enrich.num.pop_one <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop_two))
  calcLD(snp_folder=enrich_folder, popcode=2)
  enrich.pop_two <- diff.pairs
  enrich.pop_two$set <- "Enriched"
  enrich.pop_two$pop <- pop_two
  enrich.num.pop_two <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # Calculate trans-chromosomal LD stats between unenriched pathways
  cat("* Calculating trans-chromosomal SNP assocation between unenriched pathways\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop_one))
  calcLD(snp_folder=unenrich_folder, popcode=1)
  unenrich.pop_one <- diff.pairs
  unenrich.pop_one$set <- "Unenriched"
  unenrich.pop_one$pop <- pop_one
  unenrich.num.pop_one <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # population 2 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop_two))
  calcLD(snp_folder=unenrich_folder, popcode=2)
  unenrich.pop_two <- diff.pairs
  unenrich.pop_two$set <- "Unenriched"
  unenrich.pop_two$pop <- pop_two
  unenrich.num.pop_two <- diff.num
  cat(" done.\n")
  Sys.sleep(3)

  # Write out tables of SNP-SNP interactions per pathway-pathway combo
  col_pop1 <- sprintf("num_ixns_%s", pop_one)
  col_pop2 <- sprintf("num_ixns_%s", pop_two)

  ## Enriched
  ixn_enrich <- unique(enrich.pop_one[,c("pathway_pair1", "pathway_pair2")])
  enrich_df <- data.frame(interaction=ixn_enrich,
                          pop_one=enrich.num.pop_one,
                          pop_two=enrich.num.pop_two)
  colnames(enrich_df)[2] <- col_pop1
  colnames(enrich_df)[3] <- col_pop2

  write.table(enrich_df,
    sprintf("%s/enrich_num_interactions_pathway-pathway.txt", output_folder),
      col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  ixn_unenrich <- unique(unenrich.pop_one[,c("pathway_pair1", "pathway_pair2")])
  unenrich_df <- data.frame(interaction=ixn_unenrich,
                            pop_one=unenrich.num.pop_one,
                            pop_two=unenrich.num.pop_two)
  colnames(unenrich_df)[2] <- col_pop1
  colnames(unenrich_df)[3] <- col_pop2

  write.table(unenrich_df,
    sprintf("%s/unenrich_num_interactions_pathway-pathway.txt", output_folder),
      col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  # Get pathway pairs with at least 1 inter-chrom SNP-SNP association
  unenrich.num.pop_one <- unenrich.num.pop_one[which(unenrich.num.pop_one != 0)]
  unenrich.num.pop_two <- unenrich.num.pop_two[which(unenrich.num.pop_two != 0)]

  ## PLOT STATS
  cat("\n*Generating trans-chromosomal LD plots and analyses.\n")
  # Function to determine significance of inter-chrom LD SNP assocition between
  # selection-enriched pathways vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  getPvals <- function(dat, pop_name) {
    # Separate df into 'Enriched' and 'Unenriched' interactions
    enriched <- filter(dat, pop==pop_name & set=="Enriched")
    unenriched <- filter(dat, pop==pop_name & set=="Unenriched")

    for (i in 1:length(unique(enriched$ixn_num))) {
      # Subset for each enriched pathway
      enrich_path_ld[[i]] <<- subset(enriched, ixn_num==sprintf("interaction_%s", i))
      # Calculate KS pvalue for each enriched pathway against the entire
      # set of unenriched pathways; calculated per population
      ks_pvals[[i]] <<- ks.test(enrich_path_ld[[i]]$R.squared,
                                unenriched$R.squared,
                                alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Population-specific pairwise trans-chromosomal r2 per enriched pathway
  title <- "Pathway-specific SNP-SNP coevolution between enriched pathways"
  r2.xaxis.title <- "SNP-SNP association (r2)"

  # Merge all results together
  dat <- rbind(enrich.pop_one, unenrich.pop_one, enrich.pop_two, unenrich.pop_two)

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$ixn_num)))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=R.squared, colour=ixn_num, fill=ixn_num)) +
         facet_grid(set~pop) +
         geom_density(alpha=0.2) +
         xlab(r2.xaxis.title) +
         ylab(bquote(bold("Density"))) +
         scale_fill_manual(values=npal) +
         scale_colour_manual(values=npal) +
         theme(legend.position="none",
               strip.text=element_text(face="bold"))

  # 1b) eCDF plot (cumulative density at each r2)
  p2 <- ggplot(dat, aes(x=R.squared, colour=ixn_num)) +
          facet_grid(set~pop) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          scale_colour_manual(values=npal) +
          theme(legend.position="none",
                strip.text=element_text(face="bold"))

  # 3c,d) Density and eCDF at r2 > 0.2
  p3 <- p1 + xlim(0.2, 1)
  p4 <- p2 + xlim(0.2, 1)

#  cat("\n*Generating plot...")
#  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
#  title <- ggdraw() + draw_label(title, fontface="bold")
#  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
#  filename <- sprintf("%s/between_snp-snp_interchrom_ld_r2.png", output_folder)
#  save_plot(filename, both, base_height=10, base_width=13.5, base_aspect_ratio=1.2)
#  cat(sprintf(" saved to %s.\n", filename))

  # Calculate pathway-level p values per population
  enrich_path_ld <- list()
  ks_pvals <- list()
  res <- list()

  cat("* Determining significant SNP-SNP coevolution between enriched pathways.\n")
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
    cat(sprintf("** Pathway-pathway interactions with significant coevolution in %s (FDR<=0.2, N=%s):\n%s\n",
        pop, nrow(sig_paths), paste(sig_paths$pathway, collapse="\n")))
    colnames(pvals_df)[2:3] <- paste(colnames(pvals_df)[2:3], sep="_", pop)
    res[[pop]] <<- pvals_df

    # Write out results
    filename_1 <- sprintf("%s/enrich_coevolution_pval_%s.txt", output_folder, pop)
    write.table(format(pvals_df, digits=3), file=filename_1,
        col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    cat(sprintf("** Table of pathway p-values written to %s.\n", filename_1))
  }
  # Run for both populations
  invisible(mapply(calcCoev, c(pop_one, pop_two)))

  # Merge both p value dataframes
  pval_merge <- join(res[[1]], res[[2]], by=c("interaction.pathway_pair1", "interaction.pathway_pair2"))
  filename_2 <- sprintf("%s/enrich_coevolution_pval_merge_%s_%s.txt", output_folder, pop_one, pop_two)
  write.table(format(pval_merge, digits=3), file=filename_2,
      col=TRUE, row=FALSE, quote=FALSE, sep="\t")
  cat(sprintf("\n** Merged p-value tables written to %s.\n", filename_2))
}
