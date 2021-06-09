# popppaths
library(data.table)
library(plyr)
library(ggplot2)

maf <- fread("data/HM3_2010-05_phase3/HM3_pops_hg19_CEU_YRI.frq.cc", h=T, data.table=F)
fst <- fread("res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/freq/HM3_pops_hg19.fst", h=T, data.table=F)

df <- join(maf, fst)
df$dMAF <- abs((df$MAF_A)-(df$MAF_U))

pdf("../Desktop/dMAF_vs_fst.pdf")
ggplot(df, aes(dMAF, FST)) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title="dMAF vs. FST") 
dev.off()
