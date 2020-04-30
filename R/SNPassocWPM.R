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

SNPassocWPM <- function(pop_one, pop_two, snp2gene_file, ASSOC_FDR=0.05,
                        enrich_folder, unenrich_folder, output_folder) {

  # Read in snp2gene file to assign SNPs to genes
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
    ss <- c(); snp_map <- c()
    snp_assoc <- list(); pairwise_df <- list(); diff_df <- list()

    # Read PLINK files for each pathway
    if (pathway_set == "enriched")   { snp_folder <- enrich_folder }
    if (pathway_set == "unenriched") { snp_folder <- unenrich_folder }

    cat(sprintf("\nReading in %s pathways\n", pathway_set))
    bed <- list.files(path=snp_folder, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snp_folder, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snp_folder, pattern="*.fam", full.names=TRUE)

    # Skip SNP assocation calculation for pathways with <= 1 SNP
    skip <- sapply(bim, readLines)
    skip <- which(sapply(skip, length) <= 1)

    if (length(skip) > 0) {
      cat(sprintf("Skipping LD calculation for pathway %s. Not enough SNPs\n\n", skip))
      bed <- bed[-skip]
      bim <- bim[-skip]
      fam <- fam[-skip]
    } else {
      cat("All pathways have enough SNPs for LD calculation. Commencing...\n\n")
    }
    Sys.sleep(3)

    # Run SNP-SNP analysis for each pathway
    for (i in seq_along(bim)) {
      # Convert PLINK files to snpStats input format
      # Output object is a list with 3 elements ($genotypes, $fam, $map)
      # NOTE: order is important!
      path_name <- gsub("\\..*", "", basename(bim[i]))
      cat(sprintf("* Reading pathway PLINK set %s\n", path_name))
      ss[[i]] <- read.plink(bed[i], bim[i], fam[i])

      # Subset genotypes by population
      pop_code <- ifelse(population_cohort == pop_one, 1, 2)
      pop_genotypes <- which(ss[[i]]$fam$affected == pop_code)
      ss[[i]]$genotypes <- ss[[i]]$genotypes[pop_genotypes,]
      cat(sprintf("* Subsetting %i genotypes for %s population(s)\n",
        length(pop_genotypes), population_cohort))
      print(ss[[i]]$genotypes)

      # Calculate SNP association for all SNPs in pathway
      cat("* Calculating SNP-SNP association statistics\n")
      snp_assoc[[i]] <- ld(ss[[i]]$genotypes,
                           stats="R.squared",
                           depth=ncol(ss[[i]]$genotypes)-1)
      snp_map <- ss[[i]]$map # genomic location of each SNP

      # Covert association matrix into a data frame
      assoc_df <- as.matrix(snp_assoc[[i]])  # convert sparseMatrix to regular matrix
      lowerTriangle(assoc_df, byrow=FALSE, diag=TRUE) <- NA # replace lower triangle and diag with NA
      assoc_df <- melt(assoc_df)
      colnames(assoc_df)[3] <- "Rsquared"

      # Generate pariwise distance table for each SNP-SNP pair
      colnames(assoc_df)[1] <- "snp.name"
      snp_map <- subset(snp_map, select=c("snp.name", "chromosome", "position"))
      pairwise <- merge(snp_map, assoc_df, by="snp.name")
      colnames(pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
      pairwise <- merge(snp_map, pairwise, by="snp.name")
      colnames(pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2", "chr_2", "pos_2", "Rsquared")
      pairwise$dist <- abs(pairwise$pos_1 - pairwise$pos_2)
      pairwise$pathway <- path_name

      # Assign SNPs to genes
      colnames(snp2gene) <- c("snp_1", "gene_1")
      pairwise <- left_join(pairwise, snp2gene, by="snp_1")
      colnames(snp2gene) <- c("snp_2", "gene_2")
      pairwise <- left_join(pairwise, snp2gene, by="snp_2")

      # Select for SNP-SNP pairs on different chromosomes
      cat("* Subsetting for trans-chromosomal SNP-SNP pairs\n\n")
      pairwise_df[[i]] <- pairwise
      diff_df[[i]] <- filter(pairwise_df[[i]], chr_1 != chr_2)
    }

    # Generate summary files
    all_pairs_df  <- do.call("rbind", pairwise_df)
    diff_pairs_df <- do.call("rbind", diff_df)
    cat(sprintf("* Completed SNP-SNP association analysis for %i pathways\n", length(bim)))
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

  # Calculate trans-chromosomal SNP association within enriched pathways
  # per population cohort
  enrich_pop_one <- calcAssoc("enriched", pop_one)
  enrich_pop_two <- calcAssoc("enriched", pop_two)
  unenrich_pop_one <- calcAssoc("unenriched", pop_one)
  unenrich_pop_two <- calcAssoc("unenriched", pop_two)

  ## PLOT STATS
  cat("\n* Generating plots to demonstrate within-pathway SNP coevolution\n")

  # Population-specific pairwise trans-chromosomal r2 per enriched pathway
  title <- "Pathway-level SNP-SNP coevolution within enriched pathways"
  r2.xaxis.title <- "SNP-SNP association (r2)"

  # Merge all results together
  dat <- rbind(enrich_pop_one, unenrich_pop_one, enrich_pop_two, unenrich_pop_two)

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$pathway)))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=Rsquared, colour=pathway, fill=pathway)) +
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
  p2 <- ggplot(dat, aes(x=Rsquared, colour=pathway)) +
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
  plot_name <- sprintf("%s/plot_WPM_distribution_SNPassoc.png", output_folder)
  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
  title <- ggdraw() + draw_label(title, fontface="bold")
  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
  save_plot(plot_name, both, base_height=10, base_width=13.5, base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", plot_name))

  # Calculate p values per population
  enrich_path_ld <- list()
  ks_pval <- list()
  res <- list()

  cat("* Determining significant within-pathway SNP coevolution\n")

  # Determine significance of trans-chromosomal SNP-SNP association
  # per enriched pathway vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})
  calcSig <- function(population_cohort) {

    # Separate data by enriched and unenriched pathways
    enriched <- filter(dat, population == population_cohort & set == "enriched")
    unenriched <- filter(dat, population == population_cohort & set == "unenriched")

    # Prepare p-value dataframe
    pval_df <- data.frame(
      pathway=unique(enriched$pathway),
      pvalue=NA,
      qvalue=NA
    )

    # Calculate p-value per enriched pathway
    for (path in unique(enriched$pathway)) {
      test_path <- subset(enriched, pathway == path)
      ks_test <-
        suppressWarnings(
          ks.test(
            test_path$Rsquared,
            unenriched$Rsquared,
            alternative="less",
      ))
      # Insert pvalue in pval_df
      ks_pval <- ks_test$p.value
      pval_df[which(pval_df$pathway == path), "pvalue"] <- ks_pval
    }

    # FDR-correct pvalues for multiple hypothesis testing
    pval_df <- pval_df[order(pval_df$pvalue),]
    pval_df$qvalue <- p.adjust(pval_df$pvalue, method="BH")

    # Quantify number of significant pathways
    n_sig <- filter(pval_df, qvalue <= ASSOC_FDR)
    cat(sprintf("* Pathways with significant coevolution in %s population: %s\n", population_cohort, nrow(n_sig)))

    # Get summary statistics for final output dataframe
    nsnp <- enriched %>% group_by(pathway) %>% summarise(n_snps=length(Rsquared))
    stat <- enriched %>% group_by(pathway) %>% summarise(mean_r2=mean(Rsquared, na.rm=TRUE))
    perc <- enriched %>% group_by(pathway) %>% summarise(perc_99=quantile(Rsquared, 0.99, na.rm=TRUE))
    pval_df <- reduce(list(nsnp, stat, perc, pval_df), left_join, by="pathway")
    pval_df <- as.data.frame(pval_df)

    # Write out results
    fname <- sprintf("%s/table_WPM_enrich_sig_coevolution_%s.txt", output_folder, population_cohort)
    write.table(format(pval_df, digits=3), file=fname, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    cat(sprintf("* Table of pathway signifiance values written to %s\n", fname))

    # Return out results as list
    return(pval_df)
  }

  # Run for both populations
  results_pop_one <- calcSig(population_cohort = pop_one)
  results_pop_two <- calcSig(population_cohort = pop_two)
}
