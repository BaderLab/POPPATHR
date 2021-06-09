# [workdir]/ld (after running LDstats.R)
require(dplyr)
require(scales)
require(ggplot2)

dat <- readRDS("res_hc_em_groups_lc_snps_BPM/hc.diff.pairs.pop2.rds")

# check <- which(dat$chr_1 == dat$chr_2)
# length(check) #should be 0

dat <- dat[,c(1,2,7,9:13)]

# substring() removes "interaction_" string from ixn_num column to order
# the dataframe numerically by the interaction number
dat_sort <- dat[order(as.numeric(substring(dat$ixn_num, 13))),]
dat_sort <- na.omit(dat_sort)

# get dataframe listing each unique pathway x pathway interaction
pathway_ixns <- unique(dat_sort[,c(4,5)])

plot_list = list()
for (i in 1:length(unique(dat_sort$ixn_num))) {
    p = ggplot(data=subset(dat_sort, ixn_num==sprintf("interaction_%i", i)), aes(gene_1, gene_2)) +
          geom_tile(aes(fill=R.squared)) + # background colours are mapped according to the value column
          scale_fill_gradient2(low="white",
                               high=muted("midnightblue")) + # determine the colour
          theme(panel.grid.major.x=element_blank(), #no gridlines
                panel.grid.minor.x=element_blank(),
                panel.grid.major.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.background=element_rect(fill="white"), # background=white
                axis.text.x=element_blank(),
                plot.title=element_text(size=12, face="bold"),
                axis.text.y=element_blank(),
                axis.title=element_text(size=10)) +
          ggtitle(sprintf("Between-pathway interaction heatmap %i", i)) +
          theme(legend.title=element_text(face="bold", size=10)) +
          xlab(pathway_ixns[i,1]) +
          ylab(pathway_ixns[i,2]) +
          labs(fill=expression(italic(r^2)))
    plot_list[[i]] = p
}

pdf("bpm_em_group_asw.pdf")
for (i in 1:length(unique(dat_sort$ixn_num))) {
    print(plot_list[[i]])
}
dev.off()
