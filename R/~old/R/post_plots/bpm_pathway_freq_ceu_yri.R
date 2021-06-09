#[workdir]/ld/res_hc_em_groups_lc_snps_BPM_w_singletons

require(plyr)
require(ggplot2)
require(cowplot)

dat <- read.delim("between_pathway_ixns.txt", h=T)
counts <- ddply(dat, .(dat$pathway, dat$population), nrow)
names(counts) <- c("pathway", "population", "count")

count_tab <- counts[order(counts$count, decreasing=T),]

# reorder factor levels
count_tab <- transform(count_tab, pathway=factor(pathway, levels=unique(pathway)))

p <- ggplot(count_tab, aes(x=pathway, y=count, fill=population)) +
        geom_bar(stat="identity") +
        labs(x="Pathway", y="Frequency") +
        scale_fill_manual(values=c("#7D54A5", "#A6CEE3")) +
        theme(legend.position="none",
              axis.text.x=element_text(angle=45, hjust=1),
              plot.margin = unit(c(10,5,5,10),"mm"))

save_plot("bpm_pathway_freq_pop.tiff", p, base_width=10.5, base_height=6.5)
