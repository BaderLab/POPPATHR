#[workdir]/ld/res_hc_em_groups_lc_snps_BPM_w_singletons
#https://mran.microsoft.com/snapshot/2017-08-18/web/packages/ggCompNet/vignettes/examples-from-paper.html

library(network)
library(sna)
library(ggplot2)
library(dplyr)
library(tidyr)
library(geomnet)
library(scales)

## USING GGPLOT
#bpm <- read.delim("bpm_network_plot_table.txt", h=T, as.is=T)

# create plot
#set.seed(1) #keep shape consistent

#ggplot(data=bpm, aes(from_id=from_id, to_id=to_id)) +
#  geom_net(aes(colour=population),
#               layout.alg="fruchtermanreingold",
#               size=2,
#               labelon=TRUE,
#               labelgeom="text",
#               linewidth=0.8,
#               selfloops=TRUE,
#               vjust=-0.6,
#               ecolour="lightgrey",
#               directed=FALSE,
#               fontsize=3,
#               ealpha=0.5) +
#  scale_colour_manual(values = c("#7D54A5", "#A6CEE3")) +
#  xlim(c(-0.05, 1.05)) +
#  theme_net() +
#  theme(legend.position="none",
#        plot.margin = unit(c(5,20,5,20),"mm"))

## USING GGNET
library(GGally)

bpm <- read.delim("bpm_network_plot_table_num.txt", h=T, as.is=T)
bpm$from_id <- as.factor(bpm$from_id)
bpm$to_id <- as.factor(bpm$to_id)
ceu <- filter(bpm, population=="CEU")
yri <- filter(bpm, population=="YRI")

# network object
ceu_edge <- ceu[,c(1,2)]
ceu_net <- network(ceu_edge, directed=TRUE)

yri_edge <- yri[,c(1,2)]
yri_net <- network(yri_edge, directed=TRUE)

# population
#pop <- data.frame(from_id=network.vertex.names(bpm_net))
#pop <- merge(pop, bpm, by="from_id", sort=FALSE)$population
#pop <- as.factor(pop)
#bpm_net %v% "population" = as.character(pop)

# color palette
#col = c("#7D54A5", "#A6CEE3")
#names(col) = levels(pop)

p1 <- ggnet2(ceu_net,
             color="#7D54A5",
          #   alpha=0.75,
             size=8,
             label=T,
             label.size=4,
             edge.alpha=0.5,
             edge.size=0.5)

 p2 <- ggnet2(yri_net,
              color="#A6CEE3",
            #  alpha=0.75,
              size=8,
              label=T,
              label.size=4,
              edge.alpha=0.5,
              edge.size=0.5)

### BARPLOT
require(plyr)
require(cowplot)

dat <- read.delim("between_pathway_ixns.txt", h=T)
counts <- ddply(dat, .(dat$pathway, dat$population), nrow)
names(counts) <- c("pathway", "Population", "count")

count_tab <- counts[order(counts$count, decreasing=T),]

# reorder factor levels
count_tab <- transform(count_tab, pathway=factor(pathway, levels=unique(pathway)))

p3 <- ggplot(count_tab, aes(x=pathway, y=count, fill=Population)) +
        geom_bar(stat="identity") +
        labs(x="Functional pathway cluster", y="Frequency") +
        scale_fill_manual(values=c("#7D54A5", "#A6CEE3")) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              plot.margin = unit(c(10,5,5,15),"mm"))

## PLOT ALL 3 PLOTS
first = plot_grid(p1, labels="A", nrow=1, ncol=2)
second = plot_grid(p2, labels="B", nrow=1, ncol=2)
third = plot_grid(p3, labels="C")
all = plot_grid(first, second, third, labels=c('', '', ''), ncol=1)

save_plot("bpm_network_bar_freq_label2.pdf", all, base_height=14, base_width=10,
          base_aspect_ratio=1.2, ncol=1, nrow=1)
