#' Calculates SNP-SNP association statistics and perform exploratory analyses
#' within pathways
#'
#' @param pop_one (char) character code for the first population (e.g., CEU).
#' @param pop_two (char) character code for the second population (e.g., YRI).
#' @param snp2gene_file (char) path to file with snp-to-gene mappings
#'    (output of SNP2gene.R).
#' @param ASSOC_FDR (interger) FDR value to define significant pathway coevolution
#'    (default=0.05).
#' @param enrich_folder (char) path to selection-enriched pathway SNP lists.
#' @param unenrich_folder (char) path to unenriched pathway SNP lists.
#' @param output_folder (char) path to output directory.
#'
#' @return none
#' @export
#'

SNPassocBPM <- function(pop_one, pop_two, snp2gene_file, ASSOC_FDR=0.05,
                        enrich_folder, unenrich_folder, output_folder) {

  # Read in snp2gene file to assign SNPs to genes to remove any
  # matching genes in a given pathway-pathway pair
  snp2gene <- read.table(snp2gene_file, h=FALSE, as.is=TRUE)
  snp2gene <- snp2gene[,-3]

  # Calculate SNP-SNP association statistics
  # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
  # to be considered, so that depth=1 would specify calculation of LD b/w
  # immediately adjacent SNPs. For our purposes we want to determine LD
  # b/w all SNPs in each pathway despite their distance from each other,
  # so we specify depth in ld() as ((number of SNPs)-1)
  calcAssoc <- function(pathway_set, population_cohort) {

    # Initiate vectors and lists to fill up with data
    ss1 <- c(); ss2 <- c(); snp_map1 <- c(); snp_map2 <- c()
    snp_assoc <- list(); pairwise_df <- list(); diff_df <- list()

    # Read PLINK files for each pathway
    if (pathway_set == "enriched")   { snp_folder <- enrich_folder }
    if (pathway_set == "unenriched") { snp_folder <- unenrich_folder }

    cat(sprintf("\nReading in %s pathways\n", pathway_set))
    bed <- list.files(path=snp_folder, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snp_folder, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snp_folder, pattern="*.fam", full.names=TRUE)

    # Matrices of all possible pathway x pathway combinations
    bed_pair <- t(combn(bed, 2))
    bim_pair <- t(combn(bim, 2))
    fam_pair <- t(combn(fam, 2))

    # Run SNP-SNP analysis for each pathway-pathway pair
    for (i in 1:nrow(bim_pair)) {

      # Get names of each pathway in pair
      path_name1 <- gsub("\\..*", "", basename(bim_pair[i,1]))
      path_name2 <- gsub("\\..*", "", basename(bim_pair[i,2]))

      # Convert PLINK files to snpStats input format
      # Output object is a list with 3 elements ($genotypes, $fam, $map)
      # NOTE: order is important!
      cat(sprintf("BETWEEN-PATHWAY INTERACTION %i (%s total)\n", i, nrow(bim_pair)))
      cat(sprintf("* Reading SNPs in %s and %s pathways\n", path_name1, path_name2))
      ss1[[i]] <- read.plink(bed_pair[i,1], bim_pair[i,1], fam_pair[i,1])
      ss2[[i]] <- read.plink(bed_pair[i,2], bim_pair[i,2], fam_pair[i,2])

      # Subset genotypes by population
      # For first pathway of pathway x pathway interaction
      pop_code <- ifelse(population_cohort == pop_one, 1, 2)
      pop_genotypes <- which(ss1[[i]]$fam$affected == pop_code)
      ss1[[i]]$genotypes <- ss1[[i]]$genotypes[pop_genotypes,]
      cat(sprintf("\n* Subsetting %i genotypes for %s population(s)\n",
        length(pop_genotypes), population_cohort))
      print(ss1[[i]]$genotypes)

      # For second pathway of pathway x pathway interaction
      pop_code <- ifelse(population_cohort == pop_one, 1, 2)
      pop_genotypes <- which(ss2[[i]]$fam$affected == pop_code)
      ss2[[i]]$genotypes <- ss2[[i]]$genotypes[pop_genotypes,]
      cat(sprintf("\n* Subsetting %i genotypes for %s population(s)\n",
        length(pop_genotypes), population_cohort))
      print(ss2[[i]]$genotypes)

      # Calculate SNP association for all SNPs in pathway
      cat("\n* Calculating SNP-SNP association statistics\n")
      snp_assoc[[i]] <- ld(x=ss1[[i]]$genotypes, y=ss2[[i]]$genotypes, stats="R.squared")

      # Get genomic location of each SNP
      snp_map1 <- ss1[[i]]$map
      snp_map2 <- ss2[[i]]$map

      # Covert association matrix into a data frame
      assoc_df <- as.matrix(snp_assoc[[i]])  # convert sparseMatrix to regular matrix
      lowerTriangle(assoc_df, byrow=FALSE, diag=TRUE) <- NA # replace lower triangle and diag with NA
      assoc_df <- melt(assoc_df)
      colnames(assoc_df)[3] <- "Rsquared"

      # Generate pariwise distance table for each SNP-SNP pair
      colnames(assoc_df)[1:2] <- c("snp.name.1", "snp.name.2")
      snp_map1 <- subset(snp_map1, select=c("snp.name", "chromosome", "position"))
      colnames(snp_map1)[1] <- "snp.name.1"
      snp_map2 <- subset(snp_map2, select=c("snp.name", "chromosome", "position"))
      colnames(snp_map2)[1] <- "snp.name.2"
      pairwise <- merge(snp_map1, assoc_df, by="snp.name.1")
      colnames(pairwise)[1:3] <- c("snp_1", "chr_1", "pos_1")
      pairwise <- merge(snp_map2, pairwise, by="snp.name.2")
      colnames(pairwise) <- c("snp_2", "chr_2", "pos_2", "snp_1", "chr_1", "pos_1", "Rsquared")
      pairwise <- pairwise[,c(4:6, 1:3, 7)]
      pairwise$pathway_pair1 <- path_name1
      pairwise$pathway_pair2 <- path_name2
      pairwise$interaction <- sprintf("interaction_%i", i)

      # Assign SNPs to genes
      colnames(snp2gene) <- c("snp_1", "gene_1")
      pairwise <- left_join(pairwise, snp2gene, by="snp_1")
      colnames(snp2gene) <- c("snp_2", "gene_2")
      pairwise <- left_join(pairwise, snp2gene, by="snp_2")
      cat(sprintf("\tTotal SNP-SNP interactions: %i\n", nrow(pairwise)))

      # Remove any identical genes in pathway-pathway pair
      cat("* Removing any matching genes in pathway-pathway interaction... ")
      remove <- intersect(pairwise$gene_1, pairwise$gene_2)
      cat(sprintf("found %s\n", length(remove)))
      no_match <- pairwise[!(pairwise$gene_1 %in% remove),]
      no_match2 <- no_match[!(no_match$gene_2 %in% remove),]

      # Select for SNP-SNP pairs on different chromosomes
      cat("* Subsetting for trans-chromosomal SNP-SNP pairs\n")
      pairwise_df[[i]] <- no_match2
      diff_df[[i]] <- filter(pairwise_df[[i]], chr_1 != chr_2)
      cat(sprintf("\tNew number of interactions: %i\n\n", nrow(diff_df[[i]])))
    }

    # Generate summary files
    all_pairs_df  <- do.call("rbind", pairwise_df)
    diff_pairs_df <- do.call("rbind", diff_df)
    cat(sprintf("* Completed SNP-SNP association analysis for %i pathway-pathway interactions\n", nrow(bim_pair)))
    cat(sprintf("* Calculated association for %i total SNP pairs\n", nrow(all_pairs_df)))
    cat(sprintf(" --> %i total trans-chromosomal pairs\n\n", nrow(diff_pairs_df)))

    # Prepare for output
    diff_pairs_df$population <- population_cohort
    diff_pairs_df$set <- pathway_set

    # Return out
    return(diff_pairs_df)
  }

  # Calculation of all pairwise SNP combinations
  # n is the number of SNPs in each pathway, divided by 2 to account for
  # duplicate pairs ("n choose 2")
  # combos <- function(n) {
  #   return( n*(n-1)/2 ) }

  # Calculate trans-chromosomal SNP association between enriched pathways
  # per population cohort
  enrich_pop_one <- calcAssoc("enriched", pop_one)
  enrich_pop_two <- calcAssoc("enriched", pop_two)
  unenrich_pop_one <- calcAssoc("unenriched", pop_one)
  unenrich_pop_two <- calcAssoc("unenriched", pop_two)

  # Write out tables of SNP-SNP associations per pathway-pathway pair
  col_pop_one <- sprintf("n_assoc_%s", pop_one)
  col_pop_two <- sprintf("n_assoc_%s", pop_two)

  # Get number of SNP-SNP associations per pathway
  n_enrich_pop_one <- as.data.frame(table(enrich_pop_one$interaction))$Freq
  n_enrich_pop_two <- as.data.frame(table(enrich_pop_two$interaction))$Freq
  n_unenrich_pop_one <- as.data.frame(table(unenrich_pop_one$interaction))$Freq
  n_unenrich_pop_two <- as.data.frame(table(unenrich_pop_two$interaction))$Freq

  ## Enriched
  enrich_name <- sprintf("%s/table_BPM_enrich_n_snp_associations.txt", output_folder)
  enrich <- unique(enrich_pop_one[,c("pathway_pair1", "pathway_pair2")])
  enrich_df <- data.frame(enrich, pop_one=n_enrich_pop_one, pop_two=n_enrich_pop_two)
  colnames(enrich_df)[3] <- col_pop_one
  colnames(enrich_df)[4] <- col_pop_two
  write.table(enrich_df, file=enrich_name, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  ## Unenriched
  unenrich_name <- sprintf("%s/table_BPM_unenrich_n_snp_associations.txt", output_folder)
  unenrich <- unique(unenrich_pop_one[,c("pathway_pair1", "pathway_pair2")])
  unenrich_df <- data.frame(unenrich, pop_one=n_unenrich_pop_one, pop_two=n_unenrich_pop_two)
  colnames(unenrich_df)[3] <- col_pop_one
  colnames(unenrich_df)[4] <- col_pop_two
  write.table(unenrich_df, file=unenrich_name, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  ## PLOT STATS
  cat("\n* Generating plots to demonstrate between-pathway SNP coevolution\n")

  # Population-specific pairwise trans-chromosomal r2 per enriched pathway
  title <- "Pathway-level SNP-SNP coevolution between enriched pathways"
  r2.xaxis.title <- "SNP-SNP association (r2)"

  # Merge all results together
  dat <- rbind(enrich_pop_one, unenrich_pop_one, enrich_pop_two, unenrich_pop_two)

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$interaction)))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=Rsquared, colour=interaction, fill=interaction)) +
         facet_grid(set~population) +
         geom_density(alpha=0.2) +
         xlab(r2.xaxis.title) +
         ylab(bquote(bold("Density"))) +
         theme_bw() +
         scale_fill_manual(values=npal) +
         scale_colour_manual(values=npal) +
         theme(text=element_text(family="sans", size=14),
               legend.position="none",
               strip.text=element_text(face="bold"))

  # 1b) eCDF plot (cumulative density at each r2)
  p2 <- ggplot(dat, aes(x=Rsquared, colour=interaction)) +
          facet_grid(set~population) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          theme_bw() +
          scale_colour_manual(values=npal) +
          theme(text=element_text(family="sans", size=14),
                legend.position="none",
                strip.text=element_text(face="bold"))

  # 3c,d) Density and eCDF at r2 > 0.2
  p3 <- p1 + xlim(0.2, 1)
  p4 <- p2 + xlim(0.2, 1)

  cat("* Generating plot...")
  plot_name <- sprintf("%s/plot_BPM_distribution_SNPassoc.png", output_folder)
  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
  title <- ggdraw() + draw_label(title, fontface="bold")
  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
  save_plot(plot_name, both, base_height=10, base_width=13.5, base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", plot_name))

  # Calculate p values per population
  enrich_path_ld <- list()
  ks_pval <- list()
  res <- list()

  cat("* Determining significant between-pathway SNP coevolution\n")

  # Determine significance of trans-chromosomal SNP-SNP association
  # per enriched pathway vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  calcSig <- function(population_cohort) {

    # Separate data by enriched and unenriched pathways
    enriched <- filter(dat, population == population_cohort & set == "enriched")
    unenriched <- filter(dat, population == population_cohort & set == "unenriched")

    # Prepare p-value dataframe
    pval_df <- data.frame(
      interaction=unique(enriched$interaction),
      unique(enriched[,c("pathway_pair1", "pathway_pair2")]),
      pvalue=NA,
      qvalue=NA
    )

    # Calculate p-value per enriched pathway
    for (ixn in unique(enriched$interaction)) {
      test_ixn <- subset(enriched, interaction == ixn)
      ks_test <-
        suppressWarnings(
          ks.test(
            test_ixn$Rsquared,
            unenriched$Rsquared,
            alternative="less",
      ))
      # Insert pvalue in pval_df
      ks_pval <- ks_test$p.value
      pval_df[which(pval_df$interaction == ixn), "pvalue"] <- ks_pval
    }

    # FDR-correct pvalues for multiple hypothesis testing
    pval_df <- pval_df[order(pval_df$pvalue),]
    pval_df$qvalue <- p.adjust(pval_df$pvalue, method="BH")

    # Quantify number of significant pathways
    n_sig <- filter(pval_df, qvalue <= ASSOC_FDR)
    cat(sprintf("* Pathway-pathway pairs with significant coevolution in %s population: %s\n",
      population_cohort, nrow(n_sig)))

    # Write out results
    colnames(pval_df)[2:3] <- paste(colnames(pval_df)[2:3], sep="_", population_cohort)
    fname_1 <- sprintf("%s/table_BPM_enrich_sig_coevolution_%s.txt", output_folder, population_cohort)
    write.table(format(pval_df, digits=3), file=fname_1, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    cat(sprintf("* Table of pathway p-values written to %s\n", fname_1))

    # Return out results as list
    return(pval_df)
  }

  # Run for both populations
  results_pop_one <- calcSig(population_cohort = pop_one)
  results_pop_two <- calcSig(population_cohort = pop_two)

  # Merge both p value dataframes and write out
  results <- left_join(results_pop_one, results_pop_two, by="interaction")
  fname_2 <- sprintf("%s/table_BPM_enrich_sig_coevolution_%s_%s.txt", output_folder, pop_one, pop_two)
  write.table(format(results, digits=3), file=fname_2, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
  cat(sprintf("* Merged p-value table written to %s\n", fname_2))
}
