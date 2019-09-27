#' Script to determine the fraction of within and between-pathway interaction
#' discoveries supported by literature-based evidence (i.e., BioGRID)
#'
#' @param inF (char) path to BioGRID human interaction file
#'
#' @return none
#' @export
#'

######
# NOTE still need to update (2019-09-27)
######

getBiogrid <- function(inF) {
  bgrid <- read.delim(inF, h=TRUE, as.is=TRUE)

  # Filter for genetic interaction only (exclude physical interactions)
  gen_ixns <- filter(bgrid, Experimental.System.Type == "genetic")
  gen_ixns <- gen_ixns[,c(1,8,9,12:15,21,22)]

  # Replace spaces with underscore in 'Author' column to help with grep
  gen_ixns$Author <- gsub(" ", "_", gen_ixns$Author)
  outF <- sprintf("%s/hc_snps/bioGRID_interactions/BIOGRID-ORGANISM-Homo_sapiens-3.4.163_genetic_interactions.tab2", dataDir)
  write.table(gen_ixns, outF, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  # Find selection-enriched gene interactions in BioGRID
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
}



  ##
  # NOTE add nedelec eqtl mapping analysis too (easy match between genes)
  ##
