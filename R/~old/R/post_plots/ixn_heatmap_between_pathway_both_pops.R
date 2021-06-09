# [workdir]/ld (after running LDstats.R)
require(dplyr)
require(scales)
require(ggplot2)

ceu <- readRDS("res_hc_em_groups_lc_snps_BPM_w_singletons/hc.diff.pairs.pop1.rds")
yri <- readRDS("res_hc_em_groups_lc_snps_BPM_w_singletons/hc.diff.pairs.pop2.rds")

dat <- rbind(ceu, yri)
# check <- which(dat$chr_1 == dat$chr_2)
# length(check) #should be 0

# substring() removes "interaction_" string from ixn_num column to order
# the dataframe numerically by the interaction number
dat_sort <- dat[order(as.numeric(substring(dat$ixn_num, 13))),]
dat_sort <- na.omit(dat_sort)

# get dataframe listing each unique pathway x pathway interaction
pathway_ixns <- unique(dat_sort[,c(9,10)])

plot_list = list()
for (i in 1:length(unique(dat_sort$ixn_num))) {
    p = ggplot(data=subset(dat_sort, ixn_num==sprintf("interaction_%i", i)), aes(gene_1, gene_2)) +
          facet_grid(. ~ pop) +
          geom_tile(aes(fill=R.squared)) + # background colours are mapped according to the value column
          scale_fill_gradient2(low="white",
                               high=muted("midnightblue")) + # determine the colour
          theme_linedraw() +
          theme(panel.background=element_rect(fill="white"), # background=white
                panel.grid=element_blank(),
                plot.title=element_text(size=15, face="bold", hjust=0.5),
                axis.text.x=element_text(angle=90, hjust=1, vjust=1, size=6, face="bold"),
                axis.text.y=element_text(size=6, hjust=1, vjust=1, face="bold"),
                axis.title=element_text(size=12),
                axis.ticks=element_blank(),
                legend.title=element_text(face="bold", size=10)) +
          ggtitle(sprintf("%s AND %s", pathway_ixns[i,1], pathway_ixns[i,2])) +
          xlab(pathway_ixns[i,1]) +
          ylab(pathway_ixns[i,2]) +
          labs(fill=expression(italic(r^2)))
    plot_list[[i]] = p
}

pdf("bpm_em_group_ceu_yri.pdf", height=16, width=25)
for (i in 1:length(unique(dat_sort$ixn_num))) {
    print(plot_list[[i]])
}
dev.off()
