#[workdir]/ld
require(ggplot2)
require(cowplot)
require(VennDiagram)
require(gplots)

hc <- read.table("hc_snps/gseaStat_per_hc_pathway.txt", h=T)
lc <- read.table("lc_snps/gseaStat_per_lc_pathway.txt", h=T)

hc <- hc[,-c(4)]
hc <- unique(hc)

lc <- lc[,-c(4)]
lc <- unique(lc)

# get overlapping genes between hc and lc sets
same <- intersect(hc$gene, lc$gene)
# [1] 637
same <- as.data.frame(same)
colnames(same) <- "gene"

same_hc <- merge(same, hc, by="gene")
same_lc <- merge(same, lc, by="gene") # same for both
same_all <- same_hc[,-4]

# removing genes duplicated in hc and lc
hc2 <- hc[,-4]
all_hc <- rbind(hc2, same_all)
hc_unique <- all_hc[!duplicated(all_hc,fromLast = FALSE)&!duplicated(all_hc,fromLast = TRUE),]

lc2 <- lc[,-4]
all_lc <- rbind(lc2, same_all)
lc_unique <- all_lc[!duplicated(all_lc,fromLast = FALSE)&!duplicated(all_lc,fromLast = TRUE),]

same_all$set <- "Overlap"
hc_unique$set <- "Enriched"
lc_unique$set <- "Unenriched"

dat <- rbind(hc_unique, lc_unique, same_all)

# violin plot
p <- ggplot(dat, aes(x=set, y=fst, fill=set)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1, fill="white") +
        labs(y="FST", title="Pathway-level FST distributions") +
      #  stat_compare_means(label.y=1.05) +
        scale_fill_manual(values=c("#fb9a99", "gray", "#386cb0")) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank(),
              legend.position="none")

save_plot("fst_dist_enrich_unenrich_overlap_unique_snps.tiff", p,
          base_width=5, base_height=7)


# venn diagram
merge_list <- list(A = hc$gene,
                   B = lc$gene )

 venn.plot <- venn.diagram(merge_list,
                          filename="venn_enriched_unenriched_overlap.tiff",
                          fill=c("#fb9a99", "#386cb0"),
                          alpha=c(0.3,0.3),
                          main=paste("Gene overlap between the selection-enriched",
                                     "\nand unenriched pathway sets"),
                          main.cex=2.3,
                          main.fontfamily="Helvetica",
                          cat.cex=2,
                          cat.fontfamily="Helvetica",
                          cex=1.8,
                          fontfamily="Helvetica",
                          height=5000,
                          width=5000,
                          resolution=500,
                          lwd=rep(2,2),
                          scaled=TRUE
                        )

# To get the list of gene present in each Venn compartment we can use the gplots package
a <- venn(merge_list, show.plot=FALSE)

# You can inspect the contents of this object with the str() function
str(a)

# By inspecting the structure of the a object created,
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a, "intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head)
inters_df <- as.data.frame(do.call(cbind, inters))
