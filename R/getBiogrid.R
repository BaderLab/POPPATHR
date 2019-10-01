#' Script to determine the fraction of within and between-pathway interaction
#' discoveries supported by literature-based evidence (i.e., BioGRID)
#'
#' @param bgridF (char) path to BioGRID human interaction file
#' @param geneF (char) path to selection-enriched gene file (written by writePathFiles.R)
#'
#' @return none
#' @export
#'

######
# NOTE still need to update (2019-09-27)
######

getBiogrid <- function(bgridF, inDir) {

  # Reformat BioGRID annotation table
  bgrid <- read.delim(bgridF, h=TRUE, as.is=TRUE)

  # Filter for genetic interactions only (exclude physical interactions)
  gen_ixns <- filter(bgrid, Experimental.System.Type == "genetic")
  gen_ixns <- gen_ixns[,c(1,8,9,12:15,21,22)]
  gen_ixns$Author <- gsub(" ", "_", gen_ixns$Author) # replace spaces with underscore

  # Write out filtered table
  bgridF_2 <- sprintf("%s_genetic_interactions.tab2.txt", substr(inF, 0, nchar(inF)-9))
  write.table(gen_ixns, bgridF_2, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  # Find selection-enriched gene interactions in BioGRID
  outF <- sprintf("%s/sel_enriched_gen_ixns_test.txt", dirname(outF))

  # Run command
  cmd <- sprintf("grep --colour=never -wf %s %s > %s", geneF, bgridF_2, outF)
  system(cmd)

  # Annotate genes to respective selection-enriched pathway(s) to find
  # within / between-pathway interactions
  df <- read.delim(outF_2, h=F, as.is=T)
  names(df) <- names(gen_ixns)

  # analyzing sel_enriched_gen_ixns file
  dat <- read.delim("sel_enriched_gen_ixns_test.txt", h=T, as.is=T)
  table(dat$Pathway.Interaction.Type)

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
