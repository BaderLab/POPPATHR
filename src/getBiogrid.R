#' Script to determine the fraction of within and between-pathway interaction
#' discoveries supported by literature-based evidence (i.e., BioGRID)
#'
#' @param bgridF (char) path to BioGRID human interaction file.
#' @param geneF (char) path to selection-enriched gene file (written by writePathFiles.R).
#' @param outDir (char) path to write output files.
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
  bgridF_2 <- sprintf("%s_genetic_interactions.tab2.txt", substr(bgridF, 0, nchar(bgridF)-9))
  write.table(gen_ixns, bgridF_2, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

  # Find selection-enriched gene interactions in BioGRID via grep
  outF <- sprintf("%s/genes_enrich_biogrid.txt", outDir)
  cmd <- sprintf("grep --colour=never -wf %s %s > %s", geneF, bgridF_2, outF)
  system(cmd)

  # Annotate genes to respective selection-enriched pathways to find
  # within and between-pathway interactions
  df <- read.delim(outF, h=FALSE, as.is=TRUE)
  names(df) <- names(gen_ixns)


  ## STOPPED HERE
  ### Create separate table for within / between interactions
  table(dat$Pathway.Interaction.Type)
  genes <- c(dat$Official.Symbol.Interactor.A, dat$Official.Symbol.Interactor.B)
  length(unique(genes))
}

  ##
  # NOTE add nedelec eqtl mapping analysis too (easy match between genes)
  ##
