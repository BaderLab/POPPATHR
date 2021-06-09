#[workdir]/ld
require(stringr)
require(reshape2)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(ggpubr)
require(forcats) # reorder ggplot xaxis descending via fct_reorder()

# first created gseaStatFile for only hc and lc pathways
hcStatsF <- "hc_snps/gseaStatFile_hc.txt"
no_col <- max(count.fields(hcStatsF))
stats <- readLines(hcStatsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

hc_fst_list <- list()
for (i in 1:ncol(snp_stats)) {
  path <- snp_stats[-1,i]
  path <- as.data.frame(str_split_fixed(path, ",", 3))
  path[path == ""] <- NA
  path <- na.omit(path)
  hc_fst_list[i] <- as.data.frame(path[,3])
}

hc <- melt(hc_fst_list)
colnames(hc) <- c("fst", "pathway")
hc$pathway <- as.numeric(hc$pathway)
#hc$pathway <- paste("pathway", hc$pathway, sep="_")
hc$fst <- as.numeric(as.character(hc$fst))

# assign theme
hc$theme = hc$pathway
hc$theme[ hc$theme==1|hc$theme==5|hc$theme==15|hc$theme==17|hc$theme==18|
          hc$theme==19|hc$theme==20|hc$theme==22|hc$theme==23|hc$theme==26|
          hc$theme==27|hc$theme==28|hc$theme==29|hc$theme==30|hc$theme==36|
          hc$theme==37|hc$theme==38|hc$theme==46|hc$theme==49|hc$theme==51|
          hc$theme==55 ] = "Viral translation and protein targeting"
hc$theme[ hc$theme==39|hc$theme==41|hc$theme==42 ] = "Immune cell regulation"
hc$theme[ hc$theme==7|hc$theme==13|hc$theme==52 ] = "Cell chemotaxis"
hc$theme[ hc$theme==25 ] = "IL-4 mediated signaling events"
hc$theme[ hc$theme==21 ] = "Fc receptor signaling pathway"

hc$theme[ hc$theme==3|hc$theme==9|hc$theme==10|hc$theme==43|hc$theme==44|
          hc$theme==45|hc$theme==48|hc$theme==50|hc$theme==53 ] = "Growth factor and BMP signaling"
hc$theme[ hc$theme==4|hc$theme==56 ] = "Wnt/calcium signaling"
hc$theme[ hc$theme==16 ] = "Embryonic morphogenesis"
hc$theme[ hc$theme==47 ] = "Sensory organ development"

hc$theme[ hc$theme==11|hc$theme==12 ] = "Cellular response to light"
hc$theme[ hc$theme==31|hc$theme==32|hc$theme==33 ] = "Regulation of cell migration"
hc$theme[ hc$theme==2 ] = "a6b4 integrin"
hc$theme[ hc$theme==34|hc$theme==40 ] = "Regulation of lipid metabolism"
hc$theme[ hc$theme==35 ] = "Proteasomal protein catabolic process"
hc$theme[ hc$theme==6 ] = "Carbohydrate transport"

hc$theme[ hc$theme==54 ] = "Validated nuclear ER-a network"
hc$theme[ hc$theme==14 ] = "Chromatin silencing"
hc$theme[ hc$theme==24 ] = "ID signaling"
hc$theme[ hc$theme==8 ] = "Cell fate commitment"

## categorized into major themes by hand ##
hc <- read.delim("fst_pathway_em_groups.txt", h=T, as.is=T)
hc$theme <- as.factor(hc$theme)
hc$major <- as.factor(hc$major)

# make ordered factor explicit to ggplot
hc$theme <- factor(hc$theme, levels=unique(hc$theme))
hc$set <- "Selection-enriched"

# unenriched pathways for comparison (10/08/18)
lcStatsF <- "lc_snps/gseaStatFile_lc.txt"
no_col <- max(count.fields(lcStatsF))
stats <- readLines(lcStatsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

lc_fst_list <- list()
for (i in 1:ncol(snp_stats)) {
  path <- snp_stats[-1,i]
  path <- as.data.frame(str_split_fixed(path, ",", 3))
  path[path == ""] <- NA
  path <- na.omit(path)
  lc_fst_list[i] <- as.data.frame(path[,3])
}

lc <- melt(lc_fst_list)
colnames(lc) <- c("fst", "pathway")
lc$pathway <- as.numeric(lc$pathway)
#lc$pathway <- paste("pathway", lc$pathway, sep="_")
lc$fst <- as.numeric(as.character(lc$fst))
lc$theme <- "Unenriched"
lc$major <- "Unenriched"
lc$set <- "Unenriched"

both <- rbind(hc, lc)
both$set <- factor(both$set, levels=unique(both$set))

#FST boxplots per theme
#reorder boxplots by decreasing mean
p1 <- ggplot(both, aes(x=fct_reorder(theme, fst, fun=mean, .desc=TRUE), y=fst)) +
        facet_grid(.~set, scales="free_x", space="free_x") +
        geom_violin(trim=FALSE, aes(fill=set)) +
        geom_boxplot(width=0.1, fill="white") +
        guides(fill=FALSE) +
        labs(x="Functional pathway cluster",
             y=expression("F"[ST])) +
#        ggtitle("FST distribution per selection-enriched functional pathway theme") +
        scale_fill_manual(values=c("#386cb0", "#fb9a99")) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              strip.text.x=element_blank(),       #removes facet labels
              plot.margin=unit(c(5,5,5,12),"mm"))

#FST boxplots per major theme
#reorder boxplots by decreasing mean
p2 <- ggplot(both, aes(x=fct_reorder(major, fst, fun=mean, .desc=TRUE), y=fst)) +
        facet_grid(.~set, scales="free_x", space="free_x") +
        geom_violin(trim=FALSE, aes(fill=set)) +
        geom_boxplot(width=0.1, fill="white") +
        guides(fill=FALSE) +
        labs(x="Selection-enriched theme",
             y=expression(paste("Fixation index (F"[ST] ,")"))) +
#        ggtitle("FST distribution per selection-enriched functional pathway theme") +
        scale_fill_manual(values=c("#d084ae", "#dcdcdc")) +
        theme(#axis.title.y=element_blank(),
            #  axis.line.y=element_blank(),
            #  axis.ticks.y=element_blank(),
            #  axis.text.y=element_blank(),
            #  axis.text.x=element_text(angle=45, hjust=1),
              strip.text.x=element_blank(),       #removes facet labels
              plot.margin=unit(c(5,5,5,12),"mm"))

# save as pdf (2019-02-21)
save_plot("fst_compare_per_theme_plus_major_unenriched.pdf", p2,
          base_height=5.5, base_width=10)

plots <- plot_grid(p1, p2, labels=c("A", "B"), ncol=1)
save_plot("fst_per_theme_plus_major_unenriched.tiff", plots,
          base_height=11, base_width=12.5)

# add to p from `fst_violin_compare_cowplot.R`
plots <- plot_grid(p, p2, labels=c("A", "B"), nrow=1,
                   rel_widths=c(0.3,1), axis="tb", align="h")
save_plot("fst_compare_per_theme_plus_major_unenriched.png", plots,
          base_height=5.5, base_width=10)
