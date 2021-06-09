#[res]
require(tidyverse)
require(maps)
require(geosphere)
require(dplyr)
require(reshape2)
require(ggplot2)
require(cowplot)

# usa map with ceu and asw points
states <- map_data("state")
p1 <- ggplot(data=states) +
        geom_polygon(aes(x=long, y=lat, group=group), fill="#f4f4f4", color="black") +
        coord_fixed(1.3) +
        geom_point(aes(x=-110, y=34), col="#FDB462", cex=3, pch=20) + #ASW
        geom_point(aes(x=-111, y=39), col="#984EA3", cex=3, pch=20) + #CEU
        theme_nothing()

## pca
hm3 <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/pca/CEU_ASW.mds", h=T)
pnc <- read.table("out_170723_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_MAF_20-200gene/pca/CEU_ASW.mds", h=T)

hm3.pops <- read.table("../data/HM3/HM3_pops_hg19_CEU_ASW.fam", as.is=T)
hm3.pops.ctrls <- filter(hm3.pops, V6 == 1)
hm3.pops.cases <- filter(hm3.pops, V6 == 2)

pnc.pops <- read.table("../data/PNC/all/PNC_imputed_merged.CLEAN_FINAL_5sd_CEU_ASW.fam", as.is=T)
pnc.pops.ctrls <- filter(pnc.pops, V6 == 1)
pnc.pops.cases <- filter(pnc.pops, V6 == 2)

mds <- hm3
pops.ctrls <- hm3.pops.ctrls
pops.cases <- hm3.pops.cases
for (i in 1:nrow(mds)) {
 if (mds$IID[i] %in% pops.ctrls[,2]) mds[i, "population"] <- "European"
 else if (mds$IID[i] %in% pops.cases[,2]) mds[i, "population"] <- "African"
 else mds[i, "population"] <- "Other"
}
hm3.mds <- mds
hm3.mds$dataset <- "HM3"

mds <- pnc
pops.ctrls <- pnc.pops.ctrls
pops.cases <- pnc.pops.cases
for (i in 1:nrow(mds)) {
 if (mds$IID[i] %in% pops.ctrls[,2]) mds[i, "population"] <- "European"
 else if (mds$IID[i] %in% pops.cases[,2]) mds[i, "population"] <- "African"
 else mds[i, "population"] <- "Other"
}
pnc.mds <- mds
pnc.mds$dataset <- "PNC"

both <- rbind(hm3.mds[4:8], pnc.mds[4:8])
both <- na.omit(both)

## Plot results via ggplot
p2 <- ggplot(both, aes(x=C1, y=C2, colour=population)) +
       facet_wrap(~dataset) +
       geom_point(shape=4, aes(color=population)) +
       labs(x="PC1", y="PC2") +
       ylim(-0.10, 0.10) +
       scale_color_manual(values=c("#fdb462", "#984ea3", "lightgrey")) +
       scale_size_manual(values=c(2,2,0.8)) +
       theme_Publication()

plots <- plot_grid(p1, p2, labels=c("A", "B"), nrow=2)
title <- ggdraw() + draw_label("Graphical representation of European-American and African-American\npopulation structure", fontface="bold")
plots <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1))

save_plot("map_pca_ceu_asw.png", plots, base_height=7.5, base_width=6.5,
          base_aspect_ratio=1.2)
