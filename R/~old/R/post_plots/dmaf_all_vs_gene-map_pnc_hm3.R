#[res]
require(ggplot2)

hm3Dir <- "out_180215_HM3_pops_hg19_CEU-ASW_absMAF_10-300gene_500kb-dist_updated_data"
pncDir <- "out_180214_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_absMAF_CEU-ASW_10-300gene_500kb-dist"

hm3_map <- read.table(sprintf("%s/gsea/pathway_analysis/gseaStat_all_pathways.txt", hm3Dir), h=T, as.is=T)
hm3_all <- read.table(sprintf("%s/freq/markerMAF.txt", hm3Dir), h=T, as.is=T)
pnc_map <- read.table(sprintf("%s/gsea/pathway_analysis/gseaStat_all_pathways.txt", pncDir), h=T, as.is=T)
pnc_all <- read.table(sprintf("%s/freq/markerMAF.txt", pncDir), h=T, as.is=T)

hm3_map2 <- hm3_map[,c(2,3,5)]
hm3_map2$set <- "Top gene-mapped"
hm3_map2$dataset <- "HM3"
pnc_map2 <- pnc_map[,c(2,3,5)]
pnc_map2$set <- "Top gene-mapped"
pnc_map2$dataset <- "PNC"

colnames(hm3_all) <- c("snp", "dmaf")
hm3_all$set <- "All genotyped SNPs"
hm3_all$dataset <- "HM3"
colnames(pnc_all) <- c("snp", "dmaf")
pnc_all$set <- "All genotyped SNPs"
pnc_all$dataset <- "PNC"

dat <- rbind(hm3_map2, hm3_all, pnc_map2, pnc_all)

# density
p <- ggplot(dat, aes(x=dmaf)) +
        facet_grid(set~dataset, scales="free_x") +
        geom_density(colour="#7fc97f", fill="#7fc97f") +
        labs(y="Density",
             x=expression(bold(paste(Delta,"MAF")))) +
        ggtitle(expression(bold(paste("Distribution of ", Delta, "MAF statistic per gene-mapped SNP")))) +
        theme_Publication()
ggsave("dmaf_all_vs_gene_mapped.png", p, width=8, height=6)

# frequency polygraph
p <- ggplot(dat, aes(dmaf, colour=set)) +
      facet_grid(.~dataset, scales="free_y") +
      geom_freqpoly() +
      ggtitle(expression(bold(paste("Distribution of ", Delta, "MAF statistic per gene-mapped SNP")))) +
      labs(y="Number of SNPs (1000s)",
           x=expression(bold(paste(Delta,"MAF")))) +
      scale_y_continuous(labels=function(x)x/1000) + #divide scale by 1000
      scale_color_manual(values=c("#662506", "#7fc97f")) +
      theme_Publication()
ggsave("dmaf_dist_genome_top-gene-map.png", p, width=8, height=4.5)
