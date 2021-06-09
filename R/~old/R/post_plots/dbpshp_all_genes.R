# [workdir]/selection
require(splitstackshape)
require(ggplot2)

sel <- read.delim("dbPSHP_20131001.tab", h=T)
sel <- cSplit(sel, "gene", sep = ",", direction = "long")

snp2gene <- read.table("../gsea/snp2gene.txt", h=F)
genes <- as.data.frame(unique(snp2gene[,2]))
colnames(genes) <- "gene"

# selected loci HM3
gene_sel <- merge(genes, sel, by="gene")

# Remove duplicate gene entries per pathway (since one gene can have many
# positively selected loci assoicated with it)
#gene_unique <- gene_sel[!duplicated(gene_sel[1:2]),]

all_sel <- nrow(gene_unique)
unique_sel <<- length(unique(gene_unique$gene))

cat(sprintf("Total number of positively selected genes: %i.\n", all_sel))
cat(sprintf("  Number of unique genes: %i.\n", unique_sel))

hm3 <- gene_sel[,c(1,7)]
hm3$chrom <- as.character(hm3$chrom)
hm3 <- unique(hm3)
hm3_tab <- as.data.frame(table(hm3$chrom))
hm3_tab <- hm3_tab[-c(20,21),]

hm3_tab[] <- lapply(hm3_tab, function(x) if(is.factor(x)) factor(x) else x)
hm3_tab[,1] <- factor(hm3_tab[,1], levels=hm3_tab[,1][order(hm3_tab[,2])])
hm3_tab$set <- "HapMap3"

# selected loci all dbPSHP
#sel_unique <- sel[!duplicated(sel[1,3]),]
sel2 <- sel[,c(3,7)]
sel2$chrom <- as.character(sel2$chrom)
sel2 <- unique(sel2)
sel2_tab <- as.data.frame(table(sel2$chrom))
sel2_tab <- sel2_tab[-c(1,2,22,23),]

sel2_tab[] <- lapply(sel2_tab, function(x) if(is.factor(x)) factor(x) else x)
sel2_tab[,1] <- factor(sel2_tab[,1], levels=sel2_tab[,1][order(sel2_tab[,2])])
sel2_tab$set <- "dbPSHP"

both <- rbind(sel2_tab, hm3_tab)

p <- ggplot(both, aes(x=as.factor(Var1), y=Freq, fill=set)) +
        facet_grid(. ~ set, scales="free", space="free_x") +
        geom_bar(stat="identity", fill="#7fc97f") +
        labs(x="Chromosome", y="# of positively selected genes") +
        ggtitle("Distribution of positively selected loci across the humang genome") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("pos-sel_distribution_genome_unique.png", p, width=10, height=4)
