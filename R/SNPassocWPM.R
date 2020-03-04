#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#'
#' @param pop_one (char) character code for the first population.
#' @param pop_two (char) character code for the second population.
#' @param ASSOC_FDR (interger) FDR value to define significant pathway coevolution.
#' @param enrich_folder (char) path to selection-enriched pathway SNP lists
#' @param unenrich_folder (char) path to unenriched pathway SNP lists
#' @param output_folder (char) path to output directory
#'
#' @return none
#' @export
#'

SNPassocWPM <- function(pop_one, pop_two, ASSOC_FDR,
                        enrich_folder, unenrich_folder,
                        output_folder) {
  # Calculate linkage disequilbrium statistics
  # NOTE: argument 'depth' specifies the max. separation b/w pairs of SNPs
  # to be considered, so that depth=1 would specify calculation of LD b/w
  # immediately adjacent SNPs. For our purposes we want to determine LD
  # b/w all SNPs in each pathway despite their distance from each other,
  # so we specify depth in ld() as ((number of SNPs)-1)
  calcAssoc <- function(pathway_set, population) {

    # Initiate vectors and lists to fill up with data
    ss <- c(); snp_map <- c()
    snp_assoc <- list(); pairwise_df <- list(); diff_df <- list()

    # Read PLINK files for each pathway
    if (pathway_set == "enriched")   { snp_folder <- enrich_folder }
    if (pathway_set == "unenriched") { snp_folder <- unenrich_folder }

    cat(sprintf("* Reading in pathways for %s pathway set\n", pathway_set))
    bed <- list.files(path=snp_folder, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snp_folder, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snp_folder, pattern="*.fam", full.names=TRUE)

    # Skip LD calculation for pathways with <= 1 SNP
    skip <- sapply(bim, readLines)
    skip <- which(sapply(skip, length) <= 1)

    if (length(skip) > 0) {
      cat(sprintf("** Skipping LD calculation for pathway %s. Not enough SNPs\n\n", skip))
      bed <- bed[-skip]
      bim <- bim[-skip]
      fam <- fam[-skip]
    } else {
      cat("** All pathways have enough SNPs for LD calculation. Commencing...\n\n")
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
      popcode <- ifelse(population == pop_one, 1, 2)
      population <- which(ss[[i]]$fam$affected == popcode)
      population <- ifelse(popcode == 1, pop_one, pop_two)
      ss[[i]]$genotypes <- ss[[i]]$genotypes[population,]
      cat(sprintf("* Subsetting %i genotypes for %s population(s)\n", length(population), population))
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
      assoc_df <- na.omit(melt(assoc_df)) # keep all non-NA values
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

      # Select for SNP-SNP pairs on different chromosomes
      cat("* Subsetting for trans-chromosomal SNP-SNP pairs\n\n")
      pairwise_df[[i]] <- pairwise
      diff_df[[i]] <- filter(pairwise_df[[i]], chr_1 != chr_2)
    }

    # Generate summary files
    all_pairs_df  <- do.call("rbind", pairwise_df)
    diff_pairs_df <- do.call("rbind", diff_df)
    cat(sprintf("* Completed trans-chromosomal SNP association analysis for %i pathways\n", length(bim)))
    cat(sprintf("* Calculated LD for %i total SNP pairs\n", nrow(all_pairs_df)))
    cat(sprintf(" --> %i total trans-chromosomal pairs\n\n", nrow(diff_pairs_df)))

    # Prepare for output
    diff_pairs_df$population <- population
    diff_pairs_df$set <- pathway_set
    # Return out
    return(diff_pairs_df)
  }

  # Calculatin of all pairwise SNP combinations
  # n is the number of SNPs in each pathway, divided by 2 to account for
  # duplicate pairs ("n choose 2")
  # combos <- function(n) {
  #   return( n*(n-1)/2 ) }

  # Calculate trans-chromosomal SNP association within enriched pathways
  enrich_pop_one <- calcAssoc("enriched", pop_one)
  enrich_pop_two <- calcAssoc("enriched", pop_two)
  unenrich_pop_one <- calcAssoc("unenriched", pop_one)
  unenrich_pop_two <- calcAssoc("unenriched", pop_two)

  # Write out tables of SNP-SNP associations per pathway
  col_pop_one <- sprintf("n_assoc_%s", pop_one)
  col_pop_two <- sprintf("n_assoc_%s", pop_two)

  ## Get number of SNP-SNP associations per pathway
  n_enrich_pop_one <- as.data.frame(table(enrich_pop_one$pathway))$Freq
  n_enrich_pop_two <- as.data.frame(table(enrich_pop_two$pathway))$Freq
  n_unenrich_pop_one <- as.data.frame(table(unenrich_pop_one$pathway))$Freq
  n_unenrich_pop_two <- as.data.frame(table(unenrich_pop_two$pathway))$Freq

  ## Enriched
  enrich_name <- sprintf("%s/n_interactions_enriched.txt", output_folder)
  enrich <- list.files(path=enrich_folder, pattern="*.snps$", full.names=FALSE)
  enrich <- gsub("\\..*", "", enrich)
  enrich_df <- data.frame(pathway=enrich, pop_one=n_enrich_pop_one, pop_two=n_enrich_pop_two)
  colnames(enrich_df)[2] <- col_pop_one
  colnames(enrich_df)[3] <- col_pop_two
  write.table(enrich_df, file=enrich_name, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  ## Unenriched
  unenrich_name <- sprintf("%s/n_interactions_unenriched.txt", output_folder)
  unenrich <- list.files(path=unenrich_folder, pattern="*.snps$", full.names=FALSE)
  unenrich <- gsub("\\..*", "", unenrich)
  unenrich_df <- data.frame(pathway=unenrich, pop_one=n_unenrich_pop_one, pop_two=n_unenrich_pop_two)
  colnames(unenrich_df)[2] <- col_pop_one
  colnames(unenrich_df)[3] <- col_pop_two
  write.table(unenrich_df, file=unenrich_name, col=TRUE, row=FALSE, quote=FALSE, sep="\t")



  ## NOTE BEGIN EDITING HERE (2020-03-04) ######


  ## PLOT STATS
  cat("\n* Generating trans-chromosomal association plots and analyses.\n")
  # Function to determine significance of trans-chromosomal LD per each
  # enriched pathway vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})

  pval <- function(dat, population) {
    # Separate df into enriched and unenriched pathways
    enriched <- filter(dat, population==population & set=="Enriched")
    unenriched <- filter(dat, population==population & set=="Unenriched")

    for (i in 1:length(unique(enriched$pathway))) {
      # Subset for each enriched pathway
      enrich_path_ld[[i]] <<- subset(enriched, pathway==unique(enriched$pathway)[i])
      # Calculate KS pvalue for each enriched pathway against unenriched pathawys
      ks_pvals[[i]] <<- ks.test(enrich_path_ld[[i]]$Rsquared,
                               unenriched$Rsquared,
                               alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Population-specific pairwise trans-chromosomal r2 per enriched pathway
  title <- "Pathway-specific SNP-SNP coevolution within enriched pathways"
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
                strip.text = element_text(face="bold"))

  # 3c,d) Density and eCDF at r2 > 0.2
  p3 <- p1 + xlim(0.2, 1)
  p4 <- p2 + xlim(0.2, 1)

  cat("* Generating plot...")
  plot_name <- sprintf("%s/within_snp-snp_interchrom_ld_r2.png", output_folder)
  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
  title <- ggdraw() + draw_label(title, fontface="bold")
  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
  save_plot(plot_name, both, base_height=10, base_width=13.5, base_aspect_ratio=1.2)
  cat(sprintf(" saved to %s.\n", filename))

  # Calculate p values per population
  enrich_path_ld <- list()
  ks_pvals <- list()
  res <- list()

  cat("* Determining significant coevolution within enriched pathways.\n")
  calcCoev <- function(population) {
    # p value function defined above
    pval(dat=dat, population=population)

    # Generate p value dataframe
    pvals <- unlist(ks_pvals)
    pvals_df <- data.frame(pathway=enrich,
                           pval=as.data.frame(pvals),
                           fdr=p.adjust(pvals, method="BH"))

    pvals_df <- pvals_df[order(pvals_df$fdr),]
    sig_paths <- filter(pvals_df, fdr <= ASSOC_FDR)
    cat(sprintf("** Pathways with significant coevolution in %s (FDR<=%s, N=%s):\n%s\n",
        population, ASSOC_FDR, nrow(sig_paths), paste(sig_paths$pathway, collapse="\n")))
    colnames(pvals_df)[2:3] <- paste(colnames(pvals_df)[2:3], sep="_", population)
    res[[population]] <<- pvals_df

    # Write out results
    filename_1 <- sprintf("%s/enrich_coevolution_pval_%s.txt", output_folder, population)
    write.table(format(pvals_df, digits=3), file=filename_1, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
    cat(sprintf("** Table of pathway p-values written to %s.\n", filename_1))
  }
  # Run for both populations
  mapply(calcCoev, c(pop_one, pop_two))

  # Merge both p value dataframes
  pval_merge <- join(res[[1]], res[[2]], by="pathway")
  filename_2 <- sprintf("%s/enrich_coevolution_pval_merge_%s_%s.txt", output_folder, pop_one, pop_two)
  write.table(format(pval_merge, digits=3), file=filename_2, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
  cat(sprintf("\n** Merged p-value tables written to %s.\n", filename_2))
}
