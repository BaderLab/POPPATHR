#FunciSNP post analysis

### NEED TO FIX THIS SCRIPT..
### NOTE: CALCULATE FISHERS EXACT ON BAR GRAPH RES

dataDir <- paste0("/media/catherine/DATAPART1/Data/PopulationPathways/methods",
                  "/3_funciSNP/outfiles")

high_lymph_eqtl <- sprintf("%s/GTEx_eQTLs/out_2017-03-15_lymph_eQTL_high_conf", dataDir)
rand_lymph_eqtl <- sprintf("%s/GTEx_eQTLs/out_2017-03-15_lymph_eQTL_rand", dataDir)
high_blood_eqtl <- sprintf("%s/GTEx_eQTLs/out_2017-03-15_whole_blood_eQTL_high_conf", dataDir)
rand_blood_eqtl <- sprintf("%s/GTEx_eQTLs/out_2017-03-15_whole_blood_eQTL_rand", dataDir)
high_tfbs <- sprintf("%s/cistrome_TFBS_immune/out_2017-03-14_immune_TFBS_high_conf", dataDir)
rand_tfbs <- sprintf("%s/cistrome_TFBS_immune/out_2017-03-14_immune_TFBS_rand", dataDir)
high_bt_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-12_blood_Tcell_TxReg_high_conf", dataDir)
rand_bt_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-12_blood_Tcell_TxReg_rand", dataDir)
high_hb_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-13_HSC_Bcell_TxReg_high_conf", dataDir)
rand_hb_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-13_HSC_Bcell_TxReg_rand", dataDir)
high_ns_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-14_neurosph_TxReg_high_conf", dataDir)
rand_ns_txreg <- sprintf("%s/roadmap_TxReg/out_2017-03-14_neurosph_TxReg_rand", dataDir)

require(ggplot2)
require(reshape2)
################################################################################
#Get FStables for high confidence and random datasets for each test group
#3 groups: cistrome_TFBS_immune, GTEx_eQTLs, and roadmap_TxReg
res <- "FStable_all.txt"
high_lymph_eqtl <- read.table(sprintf("%s/%s", high_lymph_eqtl, res), h=T)
rand_lymph_eqtl <- read.table(sprintf("%s/%s", rand_lymph_eqtl, res), h=T)
high_blood_eqtl <- read.table(sprintf("%s/%s", high_blood_eqtl, res), h=T)
rand_blood_eqtl <- read.table(sprintf("%s/%s", rand_blood_eqtl, res), h=T)
high_tfbs <- read.table(sprintf("%s/%s", high_tfbs, res), h=T)
rand_tfbs <- read.table(sprintf("%s/%s", rand_tfbs, res), h=T)
high_bt_txreg <- read.table(sprintf("%s/%s", high_bt_txreg, res), h=T)
rand_bt_txreg <- read.table(sprintf("%s/%s", rand_bt_txreg, res), h=T)
high_hb_txreg <- read.table(sprintf("%s/%s", high_hb_txreg, res), h=T)
rand_hb_txreg <- read.table(sprintf("%s/%s", rand_hb_txreg, res), h=T)
high_ns_txreg <- read.table(sprintf("%s/%s", high_ns_txreg, res), h=T)
rand_ns_txreg <- read.table(sprintf("%s/%s", rand_ns_txreg, res), h=T)

#Filter out duplicated rows to get all unique tag SNPs
high_lymph_eqtl <- high_lymph_eqtl[!duplicated(high_lymph_eqtl$tag.snp.id), ]
rand_lymph_eqtl <- rand_lymph_eqtl[!duplicated(rand_lymph_eqtl$tag.snp.id), ]
high_blood_eqtl <- high_blood_eqtl[!duplicated(high_blood_eqtl$tag.snp.id), ]
rand_blood_eqtl <- rand_blood_eqtl[!duplicated(rand_blood_eqtl$tag.snp.id), ]
high_tfbs <- high_tfbs[!duplicated(high_tfbs$tag.snp.id), ]
rand_tfbs <- rand_tfbs[!duplicated(rand_tfbs$tag.snp.id), ]
high_bt_txreg <- high_bt_txreg[!duplicated(high_bt_txreg$tag.snp.id), ]
rand_bt_txreg <- rand_bt_txreg[!duplicated(rand_bt_txreg$tag.snp.id), ]
high_hb_txreg <- high_hb_txreg[!duplicated(high_hb_txreg$tag.snp.id), ]
rand_hb_txreg <- rand_hb_txreg[!duplicated(rand_hb_txreg$tag.snp.id), ]
high_ns_txreg <- high_ns_txreg[!duplicated(high_ns_txreg$tag.snp.id), ]
rand_ns_txreg <- rand_ns_txreg[!duplicated(rand_ns_txreg$tag.snp.id), ]

uniqueTag <- function(dat) {
    (nrow(dat)/723)*100 }

perc_eqtl <- data.frame(
      lymph=c(uniqueTag(high_lymph_eqtl), uniqueTag(rand_lymph_eqtl)),
      whole_blood=c(uniqueTag(high_blood_eqtl), uniqueTag(rand_blood_eqtl)),
      pathway_group=c("highconf", "rand"), type="eQTL"
  )

perc_tfbs <- data.frame(
      lymph=c(uniqueTag(high_tfbs), uniqueTag(rand_tfbs)),
      pathway_group=c("highconf", "rand"), type="TFBS"
  )

perc_txreg <- data.frame(
      blood_Tcell=c(uniqueTag(high_bt_txreg), uniqueTag(rand_bt_txreg)),
      HSC_Bcell=c(uniqueTag(high_hb_txreg), uniqueTag(rand_hb_txreg)),
      neurosph=c(uniqueTag(high_ns_txreg), uniqueTag(rand_ns_txreg)),
      pathway_group=c("highconf", "rand"), type="TxReg"
  )

df1 <- melt(perc_eqtl, id.vars=c("pathway_group", "type"))
df2 <- melt(perc_tfbs, id.vars=c("pathway_group", "type"))
df3 <- melt(perc_txreg, id.vars=c("pathway_group", "type"))
df.all <- rbind(df1, df2, df3)

df.all$N <- c(nrow(high_lymph_eqtl), nrow(rand_lymph_eqtl),
              nrow(high_blood_eqtl), nrow(rand_blood_eqtl),
              nrow(high_tfbs), nrow(rand_tfbs),
              nrow(high_bt_txreg), nrow(rand_bt_txreg),
              nrow(high_hb_txreg), nrow(rand_hb_txreg),
              nrow(high_ns_txreg), nrow(rand_ns_txreg))

df.all

#Plot relative percentage of tagSNPs for each bio feature category
title <- paste0("Identifying putative functional SNPs among the high-",
                "confidence\n pathway group vs. randomly selected SNPs")
ggplot(df.all, aes(x=variable, y=value, fill=pathway_group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      facet_grid(.~type, scales="free", space="free_x") +
      labs(x="Genomic features", y="Occurrence (%)") +
      ggtitle(title) +
      geom_text(aes(label=N), vjust=-0.3, size=3.5,
                    position=position_dodge(width=0.9)) +
      theme_set(theme_grey()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1))
ggsave("fs_res_all2.png", width=9, height=8)

#-------------------------------------------------------------------------------
#Function to count instances of "YES" in Promoter, Intron, Exon, and Intergenic
#colummns of each FS results table, calculates relative percentage
yesPerc <- function(dat) {
    (length(which(dat=="YES"))/length(dat))*100 }

#By genomic function
dat.high <- high_lymph_eqtl
dat.rand <- rand_lymph_eqtl

#blah <- which(dat.high$Intron == "YES" & dat.high$Intergenic == "YES")

#Create df for relative percentage of promoters, introns, exons, and intergenic
#SNPs per dataset
yes_perc <- data.frame(
      promoter=c(yesPerc(dat.high$Promoter), yesPerc(dat.rand$Promoter)),
      exon=c(yesPerc(dat.high$Exon), yesPerc(dat.rand$Exon)),
      intron=c(yesPerc(dat.high$Intron), yesPerc(dat.rand$Intron)),
      intergenic=c(yesPerc(dat.high$Intergenic), yesPerc(dat.rand$Intergenic)),
      pathway_group=c("highconf", "rand"),
      type="eQTL", tissue="lymph" )

#Reshape df for ggplot
dff1 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))
dff2 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))
dff3 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))
dff4 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))
dff5 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))
dff6 <- melt(yes_perc, id.vars=c("pathway_group", "type", "tissue"))

dff.all <- rbind(dff1, dff2, dff3, dff4, dff5, dff6)

dff.all

#Plot distribution (barplot)
title <- paste0("Identifying putative functional SNPs among the high-",
                "confidence\n pathway group vs. randomly selected SNPs")
ggplot(dff.all, aes(x=variable, y=value, fill=pathway_group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      facet_grid(~type, scales="free", space="free_x") +
    #  facet_wrap(type~tissue) +
      labs(x="Genomic function", y="Occurrence (%)") +
      ggtitle(title) +
      theme_set(theme_grey()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1))
ggsave("fs_res_func_grid.png", width=10, height=8)
