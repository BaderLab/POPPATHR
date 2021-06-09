#[res]
require(ggplot2)

dir <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/gsea/pathway_analysis/gseaStat_all_pathways.txt",
      h=T, as.is=T)
dir <- dir[,-c(4,5)]
dir <- dir[order(dir$dmaf, decreasing=T),]
colnames(dir)[3] <- paste(colnames(dir)[3], "dir", sep="_")
dir <- unique(dir)

abs <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_absMAF_20-200gene/gsea/pathway_analysis/gseaStat_all_pathways.txt",
    h=T, as.is=T)
abs <- abs[,-4]
abs <- abs[order(abs$dmaf, decreasing=T),]
colnames(abs)[3] <- paste(colnames(abs)[3], "abs", sep="_")
abs <- unique(abs)

merged <- merge(dir, abs, by=c("gene", "snp"))
merged$diff <- abs(merged$dmaf_dir - merged$dmaf_abs)

nrow(dir)
#[1] 11394
nrow(abs)
#[1] 11394
nrow(merged)
#[1] 6146
range(merged$diff)
#[1] 0.0000 0.9236
merged <- merged[order(merged$diff, decreasing=T),]

p <- ggplot(merged, aes(x=dmaf_dir, y=dmaf_abs)) +
      geom_point(alpha=0.3, colour="#7fc97f") +
      labs(x="Directional", y="Absolute") +
      ggtitle("Gene rank comparison between SNP-level\ntest statistic defintions") +
      theme_Publication()
ggsave("dmaf_directional_abs_6146-genes.png", p, width=6, height=4)
