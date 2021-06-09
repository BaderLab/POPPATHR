# Nonparametric test to compute significance of the number of interactions
# in the coevolved pathwys, relative to unenriched pathways.

# a) fraction of interactions in epistatic coevolved pathways
# b) fraction in randomly sampled pathway pairs from BIOGRID, matched for
#    size with selection-enriched pathways – do this 1000 times ---
#    sample with replacement.
# c) pvalue (fraction of times (b) is greater than (a))

# ** Make sure you adjust for the BIOGRID background. Ie. Your universe
# is all BIOGRID interactions.
# Foreground is num interactions in selection-enriched pathways.
# Background must match foreground size.

## START
# Get all pathway genes
library(fgsea)
library(reshape2)

dataDir <- "~/PopulationPathways"
source(sprintf("%s/bin/R/post_plots/readPathways.R", dataDir))

pathF <- sprintf("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt", dataDir)
path <- gmtPathways(pathF)

workDir <- sprintf("%s/res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/ld/hc_snps", dataDir)

# Remove selection-enriched pathways from pathway set
sel_path <- read.delim(sprintf("%s/results_hc.txt", workDir), h=T, as.is=T)
sel_path <- sel_path[,1]

path <- path[!names(path) %in% sel_path]

path_melt <- melt(path)
genes_unique <- as.character(unique(path_melt$value))

# Create 1000 random gene files same size as selection-enriched
sel_enrich <- read.delim(sprintf("%s/genes_unique_hc.txt", workDir), h=F, as.is=T)
permNum <- 1000

set.seed(42)
blah <- replicate(permNum, {
              sample(x=genes_unique,
                     size=nrow(sel_enrich),
                     replace=FALSE)
          })

## NOTE other method; selects all genes from single pathways instead of
## selecting genes from entire set of pathway genes
set.seed(42)
blah <- replicate(permNum, {
  i <- sample(length(path), 1)
  rand_path <- path[[i]]

  done <- FALSE
  while (!done) {
      k <- sample(length(path), 1)
      rand_path <- c(rand_path, path[[k]])
      rand_path <- unique(rand_path)

      done <- length(rand_path) >= nrow(sel_enrich)
   }

   if (length(rand_path) > nrow(sel_enrich)) {
     diff <- length(rand_path) - nrow(sel_enrich)
     j <- sample(length(rand_path), diff)

     final <- rand_path[-j]
   } else {
     final <- rand_path
   }
})

outDir <- sprintf("%s/bioGRID_interactions/permDir3", workDir)
if (!file.exists(outDir)) dir.create(outDir)

for (i in 1:ncol(blah)) {
  write.table(blah[,i],
              file=sprintf("%s/rand_path_%i.txt", outDir, i),
              col=F, row=F, quote=F)
}

ixnF <- sprintf("%s/bioGRID_interactions/BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2", workDir)

for (i in 1:permNum) {

  inF  <- sprintf("%s/rand_path_%i.txt", outDir, i)
  outF <- sprintf("%s/rand_path_match_%i.txt", outDir, i)

  if (file.exists(outF)) {
    next
  }

  cat(sprintf("file %s\n", inF))
  system(sprintf("grep -wf %s %s > %s", inF, ixnF, outF))

}

setwd("outDir")
# permDir - random sampling from pooled pathway genes (no seed set)
# permDir2 - random sampling from whole pathways (no min/max gene number)
#            removed sel-enriched pathways
# permDir3 - random sampling from pooled pathway genes (no min/max gene number)
#            removed sel-enriched pathways
#            (redundancy b/w other pathways; only 3 genes are unique
#            to the 56 sel-enriched pathways -- GREB1, LCOR, PARP14)
#            i.e., no. unique genes without removing sel-enriched = 14,298
#                  no. unique genes with removing sel-enriched =    14,295


permF <- list.files(pattern="rand_path_match*",
                    path=getwd())
perm <- lapply(permF, function(x) read.delim(x, h=F, as.is=T))
perm_ixns <- unlist(lapply(perm, nrow))

# Mean no. of mapped interactions
mean(perm_ixns)

# Get p value
enrich_ixns <- read.delim("../sel_enriched_gen_ixns_all_colour.txt", h=F, as.is=T)
perm_p <- length(which(perm_ixns > nrow(enrich_ixns)))/permNum
