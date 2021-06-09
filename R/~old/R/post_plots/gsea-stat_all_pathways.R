#[workdir]/gsea/pathway_analysis
require(stringr)
require(reshape2)
require(ggplot2)
source("../../../../bin/R/themePublication.R")

# first create easily readible gseaStat file
statsF <- "gseaStatFile.txt"
no_col <- max(count.fields(statsF))
stats <- readLines(statsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

# get 3 separate lists for gene, snp and fst (annoying but cant figure out how
# to get one loop to generate 1 list w both sets of values
#... will figure out later)
genes <- list()
snps <- list()
fst <- list()
for (i in 1:ncol(snp_stats)) {
  path <- snp_stats[-1,i]
  path <- as.data.frame(str_split_fixed(path, ",", 3))
  path[path == ""] <- NA
  path <- na.omit(path)
  genes[i] <- as.data.frame(path[,1])
  snps[i] <- as.data.frame(path[,2])
  fst[i] <- as.data.frame(path[,3])
}

all_genes <- melt(genes)[1]
colnames(all_genes) <- "gene"
all_snps <- melt(snps)[1]
colnames(all_snps) <- "snp"
all_fst <- melt(fst)
colnames(all_fst) <- c("fst", "pathway")
all_path <- cbind(all_genes, all_snps, all_fst)

all_path$fst <- as.numeric(as.character(all_path$fst))
all_path$set <- "All"
write.table(all_path, "gseaStat_all_pathways.txt", col=T, row=F, quote=F, sep="\t")

## plot fst for all snps associated with genes (GB thesis comment)
allpath <- read.table("gseaStat_all_pathways.txt", h=T, as.is=T)
hc <- read.table("../../ld/hc_snps/gseaStat_per_hc_pathway.txt", h=T, as.is=T)
lc <- read.table("../../ld/lc_snps/gseaStat_per_lc_pathway.txt", h=T, as.is=T)

dat <- rbind(allpath, hc, lc)

p <- ggplot(dat, aes(x=fst)) +
        facet_grid(. ~ set) +
        geom_density(colour="#7fc97f", fill="#7fc97f") +
        labs(y="Density", x="FST",
             title="Distribution of FST value per SNP-mapped gene") +
        theme_Publication()
ggsave("fst_all_pathways_hc_lc.png", p, width=9, height=4)
