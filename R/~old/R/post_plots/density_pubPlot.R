library(dplyr)
library(ggplot2)

dataDir <- "~/PopulationPathways/res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/ld"

## B) Within-pathway density plot
setwd(sprintf("%s/res_hc_snps_lc_snps_updated_script_R.squared", dataDir))

# pathway_7 = cell chemotaxis (sig in ceu and yri)
# pathway_16 = embryonic morphogenesis (top in ceu)
# pathway_13 = chemotaxis (top in yri)
ceu_hc <- readRDS("hc.diff.pairs.pop1.rds")
ceu_hc_path <- filter(ceu_hc, pathway == "pathway_16")
yri_hc <- readRDS("hc.diff.pairs.pop2.rds")
yri_hc_path <- filter(yri_hc, pathway == "pathway_13")
hc <- rbind(ceu_hc_path, yri_hc_path)

ceu_lc <- readRDS("lc.diff.pairs.pop1.rds")
yri_lc <- readRDS("lc.diff.pairs.pop2.rds")
lc <- rbind(ceu_lc, yri_lc)

# begin plotting
dat <- rbind(lc, hc)
dat$set <- factor(dat$set, levels=c("Unenriched", "Enriched"))

dat %>%
  group_by(set, pop) %>%
  summarize(median=median(R.squared, na.rm=TRUE))

# NOTE CEU purple = #9977ba
# NOTE YRI blue = #a6cee3

dat$colour <- NA
dat[which(dat$set=="Unenriched"),]$colour <- "#dcdcdc"
dat[which(dat$set=="Enriched" & dat$pop=="CEU"),]$colour <- "#9977ba"
dat[which(dat$set=="Enriched" & dat$pop=="YRI"),]$colour <- "#a6cee3"

p <- ggplot(dat, aes(x=R.squared, colour=set, fill=set)) +
        facet_wrap(.~pop, scales="free") +
        geom_density(alpha=0.6) +
        scale_fill_manual(values=c("#dcdcdc", "#cc79a7")) +
        scale_colour_manual(values=c("#dcdcdc", "#cc79a7")) +
      #  ggtitle("Embryonic morphogenesis (CEU) & Chemotaxis (YRI)") +
        theme_classic() +
        theme(#legend.position="none",
              text=element_text(size=20),
              axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20))

p2 <- p + xlim(0.05, 0.5)

ggsave("density_within.pdf", p, height=5.5, width=10)
ggsave("density_within_zoom.pdf", p2, height=5.5, width=10)

## D) Between-pathway density plot
setwd(sprintf("%s/res_hc_em_groups_lc_snps_BPM_w_singletons", dataDir))

# interaction = EMBRYONIC MORPHOGENESIS + IL4 MEDIATED SIGNALLING EVENTS (ceu)
# interaction = GROWTH FACTOR AND BMP SIGNALLING + REGULATION OF LIPID METABOLISM (yri)

ceu_hc <- readRDS("hc.diff.pairs.pop1.rds")
ceu_ixn <- unique(with(ceu_hc, ceu_hc[grepl("ALPHA6BETA4*", pathway_pair1) & grepl("IL4_MEDIATED*", pathway_pair2), ]$ixn_num))
ceu_hc_ixn <- filter(ceu_hc, ixn_num == ceu_ixn)

yri_hc <- readRDS("hc.diff.pairs.pop2.rds")
yri_ixn <- unique(with(yri_hc, yri_hc[grepl("GROWTH*", pathway_pair1) & grepl("LIPID*", pathway_pair2), ]$ixn_num))
yri_hc_ixn <- filter(yri_hc, ixn_num == yri_ixn)

hc <- rbind(ceu_hc_ixn, yri_hc_ixn)

ceu_lc <- readRDS("lc.diff.pairs.pop1.rds")
yri_lc <- readRDS("lc.diff.pairs.pop2.rds")
lc <- rbind(ceu_lc, yri_lc)

# begin plotting
dat <- rbind(lc, hc)
dat$set <- factor(dat$set, levels=c("Unenriched", "Enriched"))

dat %>%
  group_by(set, pop) %>%
  summarize(median=median(R.squared, na.rm=TRUE))

p <- ggplot(dat, aes(x=R.squared, colour=set, fill=set)) +
        facet_wrap(.~pop, scales="free") +
        geom_density(alpha=0.6) +
        scale_fill_manual(values=c("#dcdcdc", "#cc79a7")) +
        scale_colour_manual(values=c("#dcdcdc", "#cc79a7")) +
      #  ggtitle("Embryonic morphogenesis (CEU) & Chemotaxis (YRI)") +
        theme_classic() +
        theme(#legend.position="none",
              text=element_text(size=20),
              axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20))

p2 <- p + xlim(0.05, 0.5)

ggsave("density_between.pdf", p, height=5.5, width=10)
ggsave("density_between_zoom.pdf", p2, height=5.5, width=10)
