# [workdir]/ld/hc_snps/bioGRID_interactions
require(dplyr)

dat <- read.delim("BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt", h=T, as.is=T)

# filter for genetic interaction only (no physical interactions)
gen_ixns <- filter(dat, Experimental.System.Type == "genetic")
gen_ixns <- gen_ixns[,c(1,8,9,12:15,21,22)]

#replace spaces with underscore in 'Author' column to help with grep
gen_ixns$Author <- gsub(" ", "_", gen_ixns$Author)

write.table(gen_ixns, "BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2",
    col=T, row=F, quote=F, sep="\t")

# grep --colour=never -wf genes_unique_hc.txt BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2
# > sel_enriched_gen_ixns.txt

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
