# Script to determine pathway genes with evidence for positive selection
# via dbPSHP database (http://jjwanglab.org/dbpshp)

#' @param hcInDir (char) path to files with high-confidence pathway SNP lists
#' @param lcInDir (char) path to files with low-confidence pathway SNP lists
#' @param outDir (char) path to output directory

pathPosSel <- function(hcInDir, lcInDir, outDir) {

  # Retrieve curated evidence file
  file_url <- "ftp://jjwanglab.org/dbPSHP/curation/dbPSHP_20131001.tab"
  file_out <- "dbPSHP_20131001.tab"
  cat(sprintf("*Fetching curated evidence file from %s.\n", file_url))
  system(sprintf("curl %s -o %s/%s", file_url, outDir, file_out))

  # Fix table entries
  sel <- read.delim(sprintf("%s/%s", outDir, file_out), h=T)
  sel <- cSplit(sel, "gene", sep = ",", direction = "long")

  getSelGenes <- function(snpDir, pathSet) {

    genes <- list.files(path=snpDir, pattern="*.genes", full.names=T)
    cat("*Reading pathway gene lists...\n")
    cat(sprintf("\t%s\n", basename(genes)))

    all_genes <- list()
    for (i in 1:length(genes)) {
      all_genes[i] <- read.table(genes[i], col.names="gene") }

    gene_melt <- melt(all_genes, value.name="gene")
    colnames(gene_melt)[2] <- "pathway"
    gene_melt[,2] <- paste("pathway", gene_melt[,2], sep="_")

    all_gene <<- nrow(gene_melt)
    unique_gene <<- length(unique(gene_melt$gene))

    cat(sprintf("Total number of pathways: %i.\n", length(genes)))
    cat(sprintf("Total number of genes in pathways: %i.\n", all_gene))
    cat(sprintf("  Number of unique genes: %i.\n", unique_gene))

    cat("\n*Determining pathways with evidence for positive selection...\n")
    gene_sel <<- merge(gene_melt, sel, by="gene")

    # Remove duplicate gene entries per pathway (since one gene can have many
    # positively selected loci assoicated with it)
    # based on 'gene' and 'pubmedid' columns
    gene_unique_path <- gene_sel[!duplicated(gene_sel[1:2]),]
    gene_unique_path <- gene_unique_path[order(gene_unique_path$pathway), ]

    all_sel <- nrow(gene_unique_path)
    unique_sel <<- length(unique(gene_unique_path$gene))

    cat(sprintf("Total number of positively selected genes: %i.\n", all_sel))
    cat(sprintf("  Number of unique genes: %i.\n", unique_sel))

    outF <- sprintf("%s/pos-sel_genes_per_%s_path.txt", outDir, pathSet)
    cat(sprintf("\n*Writing evidence file to %s/%s.\n", outDir, outF))
    write.table(gene_unique_path, outF, col=T, row=F, quote=F, sep="\t")

    # Count number of positively selected genes per pathway for plotting
    sel_total <<- as.data.frame(table(gene_unique_path$pathway))
    colnames(sel_total) <<- c("pathway", "sel_genes")

    pathnames <<- basename(genes)
    pathnames <<- gsub("\\%.*", "", pathnames)
    pathnames <<- gsub("_", " ", pathnames)
  }

  sink(sprintf("%s/pathPosSel.log", outDir))

  cat(paste("\n*Determining evidence for pathway-level selection within the",
            "confidently enriched pathways...\n"))
  getSelGenes(snpDir=hcInDir, pathSet="hc")
  hc_all_num <- unique_gene
  hc_sel_num <- unique_sel
  hc_sel <- sel_total
  hc_sel$set <- "Enriched"
  hc_pathnames <- pathnames
  hc_sel[,1] <- hc_pathnames
  cat(" done.\n\n")

  cat("==================================================================\n")
  cat(paste("*Determining evidence for pathway-level selection within the",
            "nonenriched pathways...\n"))
  getSelGenes(snpDir=lcInDir, pathSet="lc")
  lc_all_num <- unique_gene
  lc_sel_num <- unique_sel
  lc_sel <- sel_total
  lc_sel$set <- "Nonenriched"
  lc_pathnames <- pathnames
  lc_sel[,1] <- lc_pathnames
  cat(" done.\n\n")

  # Combine evidence for enriched and nonenriched and plot number of
  # genes with pos selection evidence per pathway
  both <- rbind(hc_sel, lc_sel)

  # Order dataframe by increasing 'sel_genes' column for plot
  both[] <- lapply(both, function(x) if(is.factor(x)) factor(x) else x)
  both[,1] <- factor(both[,1], levels=both[,1][order(both[,2])])

  title <- "Number of positively selected genes per identified pathway"
  cat("==================================================================\n")
  cat("*Generating plot...")
  p <- ggplot(both, aes(x=pathway, y=sel_genes, fill=set)) +
          facet_grid(. ~ set, scales="free", space="free_x") +
          geom_bar(stat="identity") +
          labs(x="Pathway", y="# of positively selected genes") +
          ggtitle(title) +
          theme_Publication() +
          theme(axis.text.x = element_text(angle=45, hjust=1),
                plot.margin = unit(c(0.5,2,0.5,8),"cm")) +
          scale_fill_Publication()

  filename <- sprintf("%s/sel_genes_hc_lc.png", outDir)
  ggsave(filename, p, width=12.5, height=9)
  cat(sprintf(" saved to %s.\n", filename))

  # Significance calculation via fishers exact test
  hc_not_sel_num <- hc_all_num - hc_sel_num
  lc_not_sel_num <- lc_all_num - lc_sel_num

  cat("*Determining significance via Fisher's exact test...\n")
  mat <- matrix(c(hc_sel_num, lc_sel_num,
                 hc_not_sel_num, lc_not_sel_num),
                nrow=2,
                dimnames=list(set=c("enriched", "nonenriched"),
                              selected=c("yes", "no")))
  print(mat)
  print(fisher.test(mat, alternative="greater"))
  cat(" done.\n")

  sink()
}
