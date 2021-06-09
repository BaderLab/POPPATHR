#[workdir]/ld
require(stringr)
require(reshape2)

# first created gseaStatFile for only lc and lc pathways
hcStatsF <- "hc_snps/gseaStatFile_hc.txt"
no_col <- max(count.fields(hcStatsF))
stats <- readLines(hcStatsF)
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

hc_genes <- melt(genes)[1]
colnames(hc_genes) <- "gene"
hc_snps <- melt(snps)[1]
colnames(hc_snps) <- "snp"
hc_fst <- melt(fst)
colnames(hc_fst) <- c("fst", "pathway")
hc <- cbind(hc_genes, hc_snps, hc_fst)

hc$fst <- as.numeric(as.character(hc$fst))
hc$set <- "Enriched"

write.table(hc, "hc_snps/gseaStat_per_hc_pathway.txt", col=T, row=F, quote=F, sep="\t")

########## lc pathways
lcStatsF <- "lc_snps/gseaStatFile_lc.txt"
no_col <- max(count.fields(lcStatsF))
stats <- readLines(lcStatsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

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

lc_genes <- melt(genes)[1]
colnames(lc_genes) <- "gene"
lc_snps <- melt(snps)[1]
colnames(lc_snps) <- "snp"
lc_fst <- melt(fst)
colnames(lc_fst) <- c("fst", "pathway")
lc <- cbind(lc_genes, lc_snps, lc_fst)

lc$fst <- as.numeric(as.character(lc$fst))
lc$set <- "Nonenriched"

write.table(lc, "lc_snps/gseaStat_per_lc_pathway.txt", col=T, row=F, quote=F, sep="\t")
