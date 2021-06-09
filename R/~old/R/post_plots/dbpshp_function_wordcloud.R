# [workdir]/selection
require(splitstackshape)
require(wordcloud)

dat <- read.delim("dbPSHP_20131001.tab", h=T)
dat <- as.data.frame(cSplit(dat, "gene", sep=",", direction="long"))
gene_func <- dat[,c(3,6)]
gene_func_unique <- unique(gene_func)

funclist <- as.data.frame(table(gene_func_unique$function.))
funclist <- funclist[-1,]
funclist <- funclist[order(funclist$Freq, decreasing=T),]
funclist[,1] <- factor(funclist[,1], levels=funclist[,1][order(funclist[,2])])

nrow(funclist)
# [1] 270

# separate column by ';'
funclist <- cSplit(funclist, "Var1", sep = ";", direction = "long")
funclist <- as.data.frame(funclist)

png("dbpshp_function_wordcloud_unique2.png", width=1280, height=1280, res=300)
wordcloud(words=funclist[,1], freq=funclist[,2],
          min.freq=1, scale=c(2, 0.2),
          max.words=300, random.order=FALSE, rot.per=0.3,
          colors=brewer.pal(8, "Accent"))
dev.off()
