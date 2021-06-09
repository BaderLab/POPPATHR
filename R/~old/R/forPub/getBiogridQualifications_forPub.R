# Dropbox/manuscript/supporting_docs/
library(plyr)
library(openxlsx)
wpm <- read.delim("within_biogridID.txt", h=F)
bpm <- read.delim("between_biogridID.txt", h=F)
bio <- read.delim("~/POPPATHR/anno/BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2.txt", h=T)
names(wpm) <- names(bio)[1]
names(bpm) <- names(bio)[1]
wpm_bio <- join(wpm, bio)
bpm_bio <- join(bpm, bio)
write.xlsx(wpm_bio, "within_biogridID_expInfo.xlsx", col=T, row=F)
write.xlsx(bpm_bio, "between_biogridID_expInfo.xlsx", col=T, row=F)
