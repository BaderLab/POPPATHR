# [poppaths]/anno/baderlab
# get no. protein coding genes in anno files
require(reshape2)

# download set of protein-coding genes
pro <- read.delim("anno_info/gene_with_protein_product_20180124.txt", h=T)

# get list of baderlab anno genes
source("../../bin/R/post_scripts/readPathways.R")
# without size filter
dat <- readPathways(fname="Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                     MIN_SIZE=0, MAX_SIZE=100000000)

blah <- melt(dat)
gene_unique <- as.data.frame(unique(blah$value))
colnames(gene_unique) <- "symbol"
pro_gene <- merge(pro, gene_unique, by="symbol")

nrow(pro_gene)
# [1] 13957
nrow(gene_unique)
# [1] 14432
nrow(pro_gene) / nrow(gene_unique)
#[1] 0.967087

##############################################################################
# with size filter
dat <- readPathways(fname="Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    MIN_SIZE=20, MAX_SIZE=200)

blah <- melt(dat)
gene_unique <- as.data.frame(unique(blah$value))
colnames(gene_unique) <- "symbol"
pro_gene <- merge(pro, gene_unique, by="symbol")

nrow(pro_gene)
# [1] 13957
nrow(gene_unique)
# [1] 14432
nrow(pro_gene) / nrow(gene_unique)
