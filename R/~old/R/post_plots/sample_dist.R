require(reshape2)
require(ggplot2)

# plot with sample distributions, show how co-selection p value calculated
toy <- data.frame(v1=c(-9,-9,-9,-9,-5), v2=c(1,1,1,1,5),
                  v3=c(3,3,3,3,7), v4=c(5,5,5,5,9))

toymelt <- melt(toy, measure.vars=c("v1", "v2", "v3", "v4"))
ggplot(toymelt, aes(x=value, fill=variable)) +
      geom_density() +
      xlim(-15,15) +
      scale_fill_manual(values=c("#386cb0","#fb9a99","#fb9a99","#fb9a99")) +
      labs(x=expression(bolditalic(paste("R"^"2")))) +
      theme_classic() +
      theme(legend.position="none",
          #  axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank()) +
      theme(axis.title = element_text(face="bold"))

ggsave("toyplot_small.png")
ggsave("toyplot.png", width=4, height=1.5)

# plot of sample barplot for epistasis overlap analysis
blah <- data.frame(num=c(1,2,3,3), path=c(1,2,3,4), set="Enriched")
blah2 <- data.frame(num=c(1,1,2,3), path=c(1,2,3,4), set="Nonenriched")
dat <- rbind(blah, blah2)

ggplot(dat, aes(y=num, x=as.factor(path), fill=set)) +
    facet_grid(.~set, scales="free", space="free_x") +
    geom_bar(stat="identity", colour="black") +
    labs(x="pathway", y="gene overlap") +
    scale_fill_manual(values=c("#fb9a99","#386cb0")) +
    theme_classic() +
    theme(legend.position="none",
        #  axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme(axis.title = element_text(face="bold"))

ggsave("toyoverlap.png")

# get NES data from gsea-stat_corr_hm3_pnc.Rhistory
ggplot(nes, aes(x=value.x, y=value.y, colour=pathway_type)) +
    geom_point(aes(alpha=pathway_type)) +
    labs(x="NES (dataset 1)", y="NES (dataset 2)") +
    scale_alpha_manual(guide='none', values=list(Enriched=1, Nonenriched=1, Other=0.07)) +
    scale_color_manual(values=c("#fb9a99", "#386cb0", "lightgrey")) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    geom_vline(xintercept=0, linetype="dashed", size=0.5) +
    theme_classic() +
    theme(legend.position="none",
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank()) +
    theme(axis.title = element_text(face="bold"))

ggsave("toy_nesrep.png")
