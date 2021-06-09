# workdir
library(plyr)
library(fgsea)
library(data.table)

## SNP-FST ranks
snp_gene <- read.delim("gsea/snp2gene.txt", h=F, as.is=T)
snp_gene <- snp_gene[,-3]
colnames(snp_gene) <- c("Marker", "Gene")

snp_fst  <- read.delim("freq/markerFST.txt", h=T, as.is=T)
snp_gene_fst <- join(snp_gene, snp_fst, by="Marker")

gene_fst <- setNames(snp_gene_fst[,"CHI2"], snp_gene_fst[,"Gene"])
gene_fst <- gene_fst[order(gene_fst, decreasing=T)]

################################################################################
# Leading edge genes
lead_path_gene <- gmtPathways("gsea/pathway_analysis/gseaLEout.txt")
lead_path <- lapply(lead_path_gene, function(x) paste(x, collapse=", "))
lead_gene <- unlist(unname(lead_path))

df <- data.frame(pathway=names(lead_path), leadingEdge=lead_gene)
df$leadingEdge <- as.character(df$leadingEdge)

################################################################################
## Selection-enriched pathways
sel_enrich <- read.delim("ld/hc_snps/results_hc.txt", h=T, as.is=T)

sel_enrich_yri <- sel_enrich[,c(1:6)]
colnames(sel_enrich_yri) <- c("pathway", "size", "ES", "NES", "pval", "padj")
sel_enrich_yri <- sel_enrich_yri[order(sel_enrich_yri$NES, decreasing=T),]

## Unenriched pathways
unenrich <- read.delim("ld/lc_snps/results_lc.txt", h=T, as.is=T)

unenrich_yri <- unenrich[,c(1:6)]
colnames(unenrich_yri) <- c("pathway", "size", "ES", "NES", "pval", "padj")
unenrich_yri <- unenrich_yri[order(unenrich_yri$NES, decreasing=T),]

################################################################################
# Merge results with leading edge file
sel_lead_gene <- join(sel_enrich_yri, df, by="pathway")
sel_lead_gene <- as.data.table(sel_lead_gene)

pdf("gseaTable_selEnrich.pdf", height=15, width=25)
plotGseaTable(lead_path_gene[sel_lead_gene$pathway],
              gene_fst,
              sel_lead_gene,
              gseaParam=0.5)
dev.off()

unsel_lead_gene <- join(unenrich_yri, df, by="pathway")
unsel_lead_gene <- as.data.table(unsel_lead_gene)

pdf("gseaTable_unenrich.pdf", height=15, width=25)
plotGseaTable(lead_path_gene[unsel_lead_gene$pathway],
              gene_fst,
              unsel_lead_gene,
              gseaParam=0.5)
dev.off()

write.table(sel_lead_gene, "~/Desktop/hc_pathway_stats_leadingEdge.txt",
            col=T, row=F, quote=F, sep="\t")
write.table(unsel_lead_gene, "~/Desktop/lc_pathway_stats_leadingEdge.txt",
            col=T, row=F, quote=F, sep="\t")

################################################################################
## NOTE fgsea is not optimized for snp-level data
## since several snps can map to the same gene, causes some genes to be repeated
## in the `gene_fst` object many different times
## needs to be run at the snp level, but pathway annotations are at the gene
## level (cannot be annotated at the snp level)
pathways <- gmtPathways("../../anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt")
fgseaRes <- fgsea(pathways=pathways,
                  stats=gene_fst,
                  minSize=10,
                  maxSize=300,
                  nperm=10000)

#Warning messages:
#1: In fgsea(pathways = pathways, stats = gene_fst, minSize = 10, maxSize = 300,  :
#  There are ties in the preranked stats (54.2% of the list).
#The order of those tied genes will be arbitrary, which may produce unexpected results.
#2: In fgsea(pathways = pathways, stats = gene_fst, minSize = 10, maxSize = 300,  :
#  There are duplicate gene names, fgsea may produce unexpected results
