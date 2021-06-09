# Measures the fraction of pathways or pathway-pathway pairs with significant
# SNP coevolution had reported genetic interactions in BioGRID
# https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.4.163/BIOGRID-ORGANISM-3.4.163.tab2.zip
require(dplyr)

dataDir <- "/Users/catherineross/PopulationPathways/res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/ld"
#datF <- list.files(pattern="Homo_sapiens", path=getwd(), full.names=T) in download folder
datF <- sprintf("%s/hc_snps/bioGRID_interactions/BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt", dataDir)
dat <- read.delim(datF, h=T, as.is=T)

# Filter for genetic interaction only (no physical interactions)
gen_ixns <- filter(dat, Experimental.System.Type == "genetic")
gen_ixns <- gen_ixns[,c(1,8,9,12:15,21,22)]

# Replace spaces with underscore in 'Author' column to help with grep
gen_ixns$Author <- gsub(" ", "_", gen_ixns$Author)
outF <- sprintf("%s/hc_snps/bioGRID_interactions/BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2", dataDir)
write.table(gen_ixns, outF, col=T, row=F, quote=F, sep="\t")

# Find selection-enriched genes in BioGRID
genes <- sprintf("%s/hc_snps/genes_unique_hc.txt", dataDir)
outF_2 <- sprintf("%s/sel_enriched_gen_ixns_test.txt", dirname(outF))

# Run command
cmd <- sprintf("grep --colour=never -wf %s %s > %s", genes, outF, outF_2)
system(cmd)

# Annotate genes to respective selection-enriched pathway(s) to find
# within / between-pathway interactions
df <- read.delim(outF_2, h=F, as.is=T)
names(df) <- names(gen_ixns)



# analyzing sel_enriched_gen_ixns file
dat <- read.delim("sel_enriched_gen_ixns.txt", h=T, as.is=T)
table(dat$Pathway.Interaction.Type)
#Between        Within  Within/Between
#    71              4             41

wpm <- filter(dat, Pathway.Interaction.Type == "Within")
bpm <- filter(dat, Pathway.Interaction.Type == "Between")
wb <- filter(dat, Pathway.Interaction.Type == "Within/Between")

table(wpm$Experimental.System)
table(bpm$Experimental.System)
table(wb$Experimental.System)

genes <- c(dat$Official.Symbol.Interactor.A, dat$Official.Symbol.Interactor.B)
length(unique(genes))
#[1] 122

##
# NOTE add nedelec eqtl mapping analysis too (easy match between genes)
##
