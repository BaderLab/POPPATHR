# [workdir]/ld/stats_hc_lc_nes
# create circos plot with R
require(dplyr)
require(RCircos)

hc <- readRDS("hc.diff.pairs.rds")
toll_path <- filter(hc, pathway=="pathway_17")

toll_circos <- toll_path[,c(2,3,5,6,7)]
toll_circos$chromEnd <- toll_circos[,2]+1
toll_circos$chromEnd.1 <- toll_circos[,4]+1
colnames(toll_circos)[1:4] <- c("Chromosome", "chromStart","Chromosome.1", "chromStart.1")
toll_circos <- toll_circos[,c(1,2,6,3,4,7,5)]
toll_circos$Chromosome <- paste("chr", toll_circos$Chromosome, sep="")
toll_circos$Chromosome.1 <- paste("chr", toll_circos$Chromosome.1, sep="")

# initialize RCircos()
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- c("chrX", "chrY")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# modify and list core components
#RCircos.List.Plot.Parameters()

# initialize graphic device
out.file <- "circos_toll_pathway.pdf"
pdf(file=out.file, height=10, width=10, compress=FALSE)
RCircos.Set.Plot.Area()

# plot chromosome ideogram
RCircos.Chromosome.Ideogram.Plot()

# gene labels and connectors
# get gene data for toll-like receptor pathway
gene_loc <- read.table("../../selection/human_hg19_refGenes_01-10-2017", h=F, as.is=T)
toll_genes <- read.table("../hc_snps/TOLL-LIKE_RECEPTOR_SIGNALING_PATHWAY%GOBP%GO:0002224.genes")
colnames(toll_genes) <- "V5"
toll_genes_loc <- merge(toll_genes, gene_loc, by="V5")
toll_genes_loc <- toll_genes_loc[,c(3,4,5,1)]
colnames(toll_genes_loc) <- c("Chromosome", "chromStart", "chromEnd", "Gene")

gene.labels <- read.table("toll_pathway_gene_loc_1.txt", h=T)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(gene.labels, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(gene.labels, name.col,track.num, side)

# link plot
path.data <- toll_circos
path.data <- path.data[order(path.data$R.squared, decreasing=F),]
#lineBig <- length(which(path.data$R.squared > 0.1))

RCircos.Link.Plot(link.data=path.data, track.num=4, by.chromosome=FALSE)

# only interactions with r2 > 0.1
path.data.high <- filter(path.data, R.squared >= 0.1)
RCircos.Link.Plot(link.data=path.data.high, track.num=4, by.chromosome=FALSE)

dev.off()
