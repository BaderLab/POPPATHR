# workdir
library(plyr)
library(dplyr)
library(OneR)
library(ggplot2)

# pop-wide fst
fst <- read.delim("freq/HM3_pops_hg19.fst", h=T)
fst <- fst[,c(2,5)]
# snp-gene mapping (distance)
snpgene <- read.delim("gsea/snp2gene.txt", h=F)
colnames(snpgene) <- c("SNP", "GENE", "DIST")
# fst of snps mapped to all tested pathways (n=5998)
all_path <- read.delim("gsea/pathway_analysis/gseaStat_all_pathways.txt", h=T)
colnames(all_path)[2] <- "SNP"

# merge snp-gene + fst data
snpgene_fst <- join(snpgene, fst)
snpgene_fst <- na.omit(snpgene_fst)

# check if fst biased towards snps closer/in genes vs far away
snpgene_fst$PATHWAY <- "NO"
snpgene_fst[which(snpgene_fst$SNP %in% all_path$SNP),]$PATHWAY <- "YES"

# bin dist and get mean fst per bin
snpgene_fst$DIST_BIN <- bin(snpgene_fst$DIST, nbins=10)

# separate by in/out of pathway
yes <- filter(snpgene_fst, PATHWAY == "YES")
df_yes <- as.data.frame(with(yes, tapply(FST, DIST_BIN, mean)))
colnames(df_yes) <- "AVG_FST"
df_yes$DIST <- rownames(df_yes)
rownames(df_yes) <- NULL
df_yes$DIST <- factor(df_yes$DIST, levels=unique(df_yes$DIST))
df_yes$N <- with(yes, tapply(SNP, DIST_BIN, length)) # N per bin

no <- filter(snpgene_fst, PATHWAY == "NO")
df_no <- as.data.frame(with(no, tapply(FST, DIST_BIN, mean)))
colnames(df_no) <- "AVG_FST"
df_no$DIST <- rownames(df_no)
rownames(df_no) <- NULL
df_no$DIST <- factor(df_no$DIST, levels=unique(df_no$DIST))
df_no$N <- with(no, tapply(SNP, DIST_BIN, length)) # N per bin

#combine
df_yes$PATHWAY <- "IN_PATHWAY"
df_no$PATHWAY <- "NOT_IN_PATHWAY"
df <- rbind(df_yes, df_no)

#plot
pdf("../../../Desktop/fst_vs_snpgeneDist.pdf", width=14, height=5.5)
ggplot(df, aes(x=DIST, y=AVG_FST)) +
    facet_grid(.~PATHWAY) +
    geom_bar(stat="identity") +
    geom_text(aes(label=N), vjust=0) +
    theme_bw() +
    labs(y="FST (average)", x="SNP-gene distance", title="FST vs. SNP-gene distance") +
    theme(text=element_text(family="sans", size=14),
          axis.text.x=element_text(angle=45, hjust=1))
dev.off()
