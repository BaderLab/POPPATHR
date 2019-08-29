#' Calculate selection statistics (LD) and perform exploratory analyses
#' for two sets of variants via R snpStats package
#'
#' @param enrichDir (char) path to selection-enriched pathway SNP lists
#' @param unenrichDir (char) path to unenriched pathway SNP lists
#' @param pop1 (char) character code for the first population (controls).
#' @param pop2 (char) character code for the second population (cases).
#' @param outDir (char) path to output directory
#'
#' @return none
#' @export
#'

LDstatsWPM <- function(enrichDir, unenrichDir, pop1, pop2, outDir) {
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
    bed <- list.files(path=snpDir, pattern="*.bed", full.names=TRUE)
    bim <- list.files(path=snpDir, pattern="*.bim", full.names=TRUE)
    fam <- list.files(path=snpDir, pattern="*.fam", full.names=TRUE)

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
    Sys.sleep(3)

    # Convert PLINK files to snpStats input format
    # Output object is a list with 3 elements ($genotypes, $fam, $map)
    # NOTE: order is important!
    for (i in seq_along(bim)) {
      path.name <- substr(basename(bim[i]), 0, nchar(basename(bim[i]))-4)
      cat(sprintf("*Reading pathway %s PLINK set\n", path.name))

      ss[[i]] <- read.plink(bed[i], bim[i], fam[i])

      # Subset genotypes by population
      pop <- which(ss[[i]]$fam$affected == popcode)
      ss[[i]]$genotypes <- ss[[i]]$genotypes[pop, ]

      cat(sprintf("*Subsetting %i genotypes based on %s population(s)\n",
          length(pop),
          if (popcode == 1) {pop1}
          else if (popcode == 2) {pop2} ))
      print(ss[[i]]$genotypes)

      cat("*Calculating LD statistics\n")
      ld.calc[[i]] <- ld(ss[[i]]$genotypes,
                         stats="R.squared",
                         depth=ncol(ss[[i]]$genotypes)-1)
      snp.map <- ss[[i]]$map # genomic location of each SNP

      # Turn each LD matrix into a data frame
      r2 <- as.matrix(ld.calc[[i]])  # convert sparseMatrix to regular matrix
      lowerTriangle(r2, byrow=FALSE, diag=TRUE) <- NA # replace lower triangle and diag with NA
      r2 <- na.omit(melt(r2)) # keep all non-NA values
      colnames(r2)[3] <- "R.squared"

      # Generate pariwise distance table for each SNP-SNP pair
      colnames(r2)[1] <- "snp.name"
      snp.map <- subset(snp.map, select=c("snp.name", "chromosome", "position"))

      pairwise <- merge(snp.map, r2, by="snp.name")
      colnames(pairwise)[1:4] <- c("snp_1", "chr_1", "pos_1", "snp.name")
      pairwise <- merge(snp.map, pairwise, by="snp.name")
      colnames(pairwise) <- c("snp_1", "chr_1", "pos_1", "snp_2",
                              "chr_2", "pos_2", "R.squared")
      pairwise$dist <- abs(pairwise$pos_1 - pairwise$pos_2)
      pairwise$pathway <- path.name

      # NOTE: rounding may affect p value caluclation by KS test (03/08/2018)
      #pairwise.df[[i]] <- pairwise %>% mutate(R.squared=round(R.squared, 3))
      #pairwise.df[[i]] <- pairwise %>% mutate(D.prime=round(D.prime, 3))
      pairwise.df[[i]] <- pairwise %>% mutate(R.squared=R.squared)

      cat("*Subsetting for interchromosomal SNP-SNP pairs\n")
      diff[[i]] <- filter(pairwise.df[[i]], chr_1 != chr_2)
      diff.r2[[i]] <- select(diff[[i]], R.squared) %>% unlist
      cat("done.\n\n")
    }

    all.pairs <<- do.call("rbind", pairwise.df)
    diff.pairs <<- do.call("rbind", diff)
    diff.r2.mean <<- sapply(diff.r2, mean) # mean r2 value per pathway
    diff.num <<- sapply(diff.r2, length)   # all SNP-SNP pairs per path

    cat(sprintf("*Finished inter-chr LD analysis for %i pathways.\n", length(bim)))
    cat(sprintf("*Calculated LD for %i total SNP pairs.\n", nrow(all.pairs)))
    cat(sprintf(" --> %i total interchromosomal pairs.\n\n", sum(diff.num)))
  }

  # How all pairwise combinations are calculated
  # n is the number of SNPs in each pathway, divided by 2 to account for
  # duplicate pairs ("n choose 2")
  # combos <- function(n) {
  #   return( n*(n-1)/2 ) }

  # Calculate inter-chromosomal LD stats within enriched pathways
  cat("=======================================================================")
  cat("**Measuring trans-chromosomal LD within selection-enriched pathways\n")

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

  # Calculate inter-chromosomal LD stats within unenriched pathways
  cat("=======================================================================")
  cat("**Measuring trans-chromosomal LD within unenriched pathways\n")

  # population 1 genotypes only
  cat(sprintf("\n========== AMONG %s GENOTYPES ==========\n", pop1))
  calcLD(snpDir=unenrichDir, popcode=1)
  unenrich.pop1 <- diff.pairs
  unenrich.pop1$set <- "Unenriched"
  unenrich.pop1$pop <- pop1
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

  # Write out tables of SNP-SNP interactions per pathway
  col_pop1 <- sprintf("num_ixns_%s", pop1)
  col_pop2 <- sprintf("num_ixns_%s", pop2)

  ## Enriched
  path_enrich <- list.files(path=enrichDir, pattern="*.snps$", full.names=FALSE)
  path_enrich <- gsub("\\%.*", "", path_enrich)
  path_enrich <- substr(path_enrich, 0, nchar(path_enrich)-5)
  enrich_df <- data.frame(pathway=path_enrich,
                          pop1=enrich.num.pop1,
                          pop2=enrich.num.pop2)
  colnames(enrich_df)[2] <- col_pop1
  colnames(enrich_df)[3] <- col_pop2

  write.table(enrich_df, sprintf("%s/enrich_num_interactions.txt", outDir),
      col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  # Unenriched
  path_unenrich <- list.files(path=unenrichDir, pattern="*.snps$", full.names=FALSE)
  path_unenrich <- gsub("\\%.*", "", path_unenrich)
  path_unenrich <- substr(path_unenrich, 0, nchar(path_unenrich)-5)
  unenrich_df <- data.frame(pathway=path_unenrich,
                            pop1=unenrich.num.pop1,
                            pop2=unenrich.num.pop2)
  colnames(unenrich_df)[2] <- col_pop1
  colnames(unenrich_df)[3] <- col_pop2

  write.table(unenrich_df, sprintf("%s/unenrich_num_interactions.txt", outDir),
      col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  ## PLOT STATS
  cat("\n*Generating inter-chromosomal LD plots and analyses.\n")
  # Function to determine significance of inter-chromosomal LD per each
  # enriched pathway vs. cumulative set of unenriched pathways
  # via the KS test (alternative=less{the CDF of x lies below+right of y})

  getPvals <- function(dat, pop_name) {
    # Separate df into 'enriched' and 'unenriched' pathways
    enriched <- filter(dat, pop==pop_name & set=="Enriched")
    unenriched <- filter(dat, pop==pop_name & set=="Unenriched")

    for (i in 1:length(unique(enriched$pathway))) {
      # Subset for each enriched pathway
      enrich_path_ld[[i]] <<- subset(enriched, pathway==unique(enriched$pathway)[i])
      # Calculate KS pvalue for each enriched pathway against unenriched pathawys
      ks_pvals[[i]] <<- ks.test(enrich_path_ld[[i]]$R.squared,
                               unenriched$R.squared,
                               alternative="less")
      ks_pvals[[i]] <<- ks_pvals[[i]]$p.value
    }
  }

  # Population-specific pairwise inter-chromosomal r2 per enriched pathway
  title <- "Pathway-specific SNP-SNP coevolution within enriched pathways"
  r2.xaxis.title <- "LD value per SNP-SNP pair"

  # Merge all results together
  dat <- rbind(enrich.pop1, unenrich.pop1, enrich.pop2, unenrich.pop2)

  # Set colour palette and number of colours needed
  cols <- colorRampPalette(brewer.pal(8, "Accent"))
  npal <- cols(length(unique(dat$pathway)))

  # 1a) Density distribution plot
  p1 <- ggplot(dat, aes(x=R.squared, colour=pathway, fill=pathway)) +
         facet_grid(set ~ pop) +
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
  p2 <- ggplot(dat, aes(x=R.squared, colour=pathway)) +
          facet_grid(set ~ pop) +
          stat_ecdf() +
          xlab(r2.xaxis.title) +
          ylab(bquote(bold("Cumulative density"))) +
          theme_bw() +
          scale_colour_manual(values=npal) +
          theme(text=element_text(family="sans", size=14),
                legend.position="none",
                strip.text = element_text(face="bold"))

  # 3c,d) Density and eCDF at x-axis > 0.2
  p3 <- p1 + xlim(0.2, 1)
  p4 <- p2 + xlim(0.2, 1)

  cat("*Generating plot...")
  both <- plot_grid(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
  title <- ggdraw() + draw_label(title, fontface="bold")
  both  <- plot_grid(title, both, ncol=1, rel_heights=c(0.1, 1))
  filename <- sprintf("%s/within_snp-snp_interchrom_ld_r2.png", outDir)
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
    pvals_df <- data.frame(pathway=path_enrich,
                           pval=as.data.frame(pvals),
                           fdr=p.adjust(pvals, method="BH"))

    pvals_df <- pvals_df[order(pvals_df$fdr),]
    sig_paths <- filter(pvals_df, fdr <= 0.2)
    cat(sprintf("Pathways with significant coevolution in %s (FDR<=0.2, N=%s):\n%s\n",
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
  pval_merge <- join(res[[1]], res[[2]], by="pathway")
  filename_2 <- sprintf("%s/enrich_coevolution_pval_merge_%s_%s.txt", outDir, pop1, pop2)
  write.table(format(pval_merge, digits=3), file=filename_2,
      col=TRUE, row=FALSE, quote=FALSE, sep="\t")
  cat(sprintf("\n*Merged p-value tables written to %s.\n", filename_2))
}
