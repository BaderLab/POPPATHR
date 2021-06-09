#[poppaths]
# Filter out genes from pathway .gmt file (e.g., only keep those w/evidence for
# positive selection)
require(tools)
source("bin/R/post_plots/readPathways.R")

annoDir <- "anno/baderlab"
selDir  <- "res/out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/selection"

gmtFile <- sprintf("%s/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
              annoDir)
myGenes <- read.table(sprintf("%s/dbPSHP_20131001_genes_only.txt", selDir),
              h=F, as.is=T)
myGenes <- as.character(unlist(myGenes))

pathList <- readPathways(gmtFile, MIN_SIZE=0, MAX_SIZE=100000)
for (k in names(pathList)) {
   pathList[[k]] <- intersect(pathList[[k]], myGenes);
  if (length(pathList[[k]]) < 1) pathList[[k]] <- NULL
}

# write file
outF <- sprintf("%s_selection_genes.gmt", file_path_sans_ext(gmtFile))
sink(outF)
for (k in names(pathList)) {
  cat(sprintf("%s\t%s\t%s\n", k, k, paste(pathList[[k]],collapse="\t")))
}
sink()
