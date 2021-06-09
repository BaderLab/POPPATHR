# Format CORUM protein complex data for GSEA pathway annotation input

file <- "http://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt"
annoDir <- "/media/catherine/DATAPART1/Data/PopulationPathways/anno/CORUM"

system(sprintf("wget -P %s %s", annoDir, file))
dat <- read.delim("coreComplexes.txt", h=T, as.is=T)
human <- which(dat$Organism == "human") #filter for 'Human'
human_dat <- dat[human, ]
human_dat$IDnew <- paste(human_dat$ComplexID, human_dat$ComplexName, sep="_")
human_dat <- human_dat[,c(20, 10, 13)] #IDnew, GO.description, subunits.Gene.name.
library(splitstackshape)
protein_complex <- cSplit(human_dat, splitCols = "subunits.Gene.name.",
                          sep = ";", direction = "wide", drop = FALSE) #split genes into separate columns
protein_complex <- protein_complex[,-3]
protein_complex <- sapply(protein_complex, as.character)
protein_complex[is.na(protein_complex)] <- ""
protein_complexdf <- as.data.frame(protein_complex)
write.table(protein_complexdf, file=sprintf("%s/coreComplexes_pathwayanno.tab",
            annoDir), col=F, row=F, quote=F, sep="\t")
