require(data.table)
require(reshape2)
require(ggplot2)
require(stringr)
require(tidyr)
require(dplyr)

allperm <- list.files(pattern="*.txt$", full.names=F)
permsample <- sample(allperm, 10)

permlist <- list()
for (i in 1:length(permsample)) {
 permlist[[i]] <- fread(permsample[i], h=T, data.table=F)
 permlist[[i]]$perm <- sprintf("perm_%i", i)
}

realperm <- fread("../gsea/freq/markerMAF.txt", h=T, data.table=F)
realperm$perm <- "real"
permlist[[11]] <- realperm

pdf("perm_maf.pdf")
d <- density(as.numeric(as.character(permlist[[11]]$CHI2)))
plot(d, main=sprintf("Delta MAF - real", i))

for (i in 1:10) {
  d <- density(as.numeric(as.character(permlist[[i]]$CHI2)))
  plot(d, main=sprintf("Delta MAF - random permutation %i", i))
}
dev.off()

################################################################################
## pathway-level perm plots
hcStatsF <- "gseaStatFile.txt"
no_col <- max(count.fields(hcStatsF))

hcStats <- read.table(hcStatsF, sep="\t", fill=T, col.names=1:no_col)
#random_paths <- hcStats[sample(nrow(hcStats), 6), ] #get 6 random pathways
stats <- t(hcStats) #transpose data

permstats <- fread("markerMAF_perm.txt", h=T, data.table=F)
charvec <- sprintf("perm.%i", seq(1:100))

pathlist <- list()
for (i in 1:ncol(stats)) {
 pathlist[[i]] <- stats[-1,i]
 pathlist[[i]] <- as.data.frame(str_split_fixed(pathlist[[i]], ",", 3))
 pathlist[[i]][pathlist[[i]] == ""] <- NA
 pathlist[[i]] <- na.omit(pathlist[[i]])
 pathlist[[i]]$path <- sprintf("path_%i", i)
 colnames(pathlist[[i]]) <- c("Gene", "Marker", "real.stat", "path")

 #add perm stats (split into 100 columns)
 pathlist[[i]] <- merge(pathlist[[i]], permstats, by="Marker")
 pathlist[[i]] <- pathlist[[i]] %>% separate(CHI2_PERM, charvec, ",")
}

path.df <- do.call("rbind", pathlist)
path.df <- path.df[,c(1,2,4,3,5:ncol(path.df))]
saveRDS(path.df, "pathway_mafstats_real_perm.rds")

for (i in 4:ncol(path.df)) {
  path.df[,i] <- as.numeric(as.character(path.df[,i])) #convert values to numeric
}

pdf("pathway_boxplots.pdf", width=19, height=8)
path_num <- length(unique(path.df$path))
for (i in 1:path_num) {
  blah <- filter(path.df, path==sprintf("path_%i", i))
  boxplot(blah[,4:ncol(blah)], main=sprintf("Pathway %i", i),
          ylab=expression(paste(Delta, "MAF")), las=2)
  abline(h=0, col="red")
}
dev.off()

################################################################################
## pathway-level maf distribution per snp
ggplot(maf.high.final, aes(x=SNP, y=MAF, colour=POP, fill=POP)) +
  geom_point() +
  ggtitle("MAF distribution - DNA CONFORMATION CHANGE (1KGP)") +
  geom_line(aes(group=SNP), colour="grey") +
  theme_set(theme_minimal()) +
  theme(plot.title=element_text(hjust=0.5),
        legend.position="top",
        legend.title=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1))
