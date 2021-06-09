#!/usr/bin/Rscript

# Removes top GSEA genes from a snp to gene map file, determined via pathways
# with an FDR < 0.3

#Input files
dataset <- "HM3"
dataDir <- "/media/catherine/DATAPART1/Data/PatientNetworks/Catherine"
inputDir  <- sprintf("%s/GWAS_GSEA/%s/CEU_ASW/out_160929_MAFdiff_CEU-ASW",
        dataDir, dataset)
snp2geneDir <- sprintf("%s/GWAS_GSEA/%s/snp2gene", dataDir, dataset)
dt	  <- format(Sys.Date(),"%y%m%d")


GSEAstatF  <- sprintf("%s/gseaStatFile.txt", inputDir)
GSEAresF <- sprintf("%s/results.txt", inputDir)
snp2geneF <- sprintf("%s/snp2gene.txt", snp2geneDir)
#logF  <- file(sprintf("%s/removeTopGenes.log", snp2geneDir), open="a")
#------------------------------------------------------------------------------
removeTopGenes <- function(snp2geneF, GSEAresF, GSEAstatF,
    snp2geneDir){

  # Read in snp2gene file
  allGenes  <- read.delim(snp2geneF, h=F, as.is=T, sep="\t")
  colnames(allGenes)  <- c("snp", "gene", "bp")
  print("Reading snp2gene file...", file=logF, sep="\n")
  print(rapply(allGenes, function(x)length(unique(x))), file=logF, sep="\n")

  # Read in GSEA results and GSEA statistics file
  GSEAres   <- read.delim(GSEAresF, h=T, as.is=T, sep="\t")
  GSEAstat  <- read.delim(GSEAstatF, h=F, as.is=T, sep="\t")
  GSEAstat  <- subset(GSEAstat, select=c(V1, V2))
  colnames(GSEAstat) <- c("Geneset", "stats")

  # Determine top genes from GSEA results file (FDR < 0.3)
  top_GSEAres <- subset(GSEAres, GSEAres$FDR < 0.3)

  # Merge top genes with gseaStatFile
  top_GSEAres_geneset <- subset(top_GSEAres, select=Geneset)
  top_GSEAstat <- merge(top_GSEAres_geneset, GSEAstat, by="Geneset")

  # Separate stat column of 'GSEAstat' into 3 columns
  library(stringr)
  top_GSEAstat <- as.data.frame(str_split_fixed(
          as.character(top_GSEAstat$stats),",", 3))
  colnames(top_GSEAstat) <- c("gene", "snp", "MAFdiff")

  # Remove all genes in snp2gene from top_GSEAstat
  ##genes_to_remove <- which(top_GSEAstat$gene %in% allGenes$gene)
  genes_to_remove <- which(allGenes$gene %in% top_GSEAstat$gene)
  genes_to_remove <- (allGenes[genes_to_remove,])

## SP added: one line to remove top gsea genes from allgenes
## this line updates allGenes
## allGenes <- allGenes[-which(allGenes$gene %in% top_GSEAstat$gene),]
## this line gets the genes that are in top_GSEAstat and in allGenes
## gone_genes <- allGenes$gene[which(allGenes$gene %in% top_GSEAstat$gene)]
## now print the unique list of these genes, alphabetically sorted so you can spot check
## sort(unique(gone_genes))
## cat(sprintf("Removed genes: %s", paste(gone_genes, collapse=",")))

    colnames(genes_to_remove)  <- c("snp", "gene", "bp")
    print("Determining genes to remove from snp2gene file (pathway FDR < 0.3)...")
    print(rapply(genes_to_remove, function(x)length(unique(x))))
    write.table(genes_to_remove, file=sprintf("%s/snp2gene_removedGenes.txt",
          snp2geneDir), col=F, row=F, sep="\t", quote=F)

    removed_genes <- allGenes[ ! allGenes$gene %in% genes_to_remove$gene, ]
    colnames(removed_genes)  <- c("snp", "gene", "bp")
    print("Total number of genes in new snp2gene file..")
    print(rapply(removed_genes, function(x)length(unique(x))))

    print("Writing out new snp2gene file..")
    write.table(removed_genes, file=sprintf("%s/snp2gene_noTopGenes.txt",
          snp2geneDir), col=F, row=F, sep="\t", quote=F)
    print("Done.")
}

removeTopGenes(snp2geneF, GSEAresF, GSEAstatF, snp2geneDir)
