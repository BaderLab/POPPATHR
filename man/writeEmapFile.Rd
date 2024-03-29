% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeEmapFile.R
\name{writeEmapFile}
\alias{writeEmapFile}
\title{Writes enrichment file for use in Cytoscape EnrichmentMap}
\usage{
writeEmapFile(results_file, ENRICH_NES = 3, out_file)
}
\arguments{
\item{results_file}{(list) file paths to results.txt for each population comparison
(generated by setupGSEArun.R).}

\item{ENRICH_NES}{(integer) NES cutoff to select validated selection-enriched
pathways (default=3).}

\item{out_file}{(char) name of output file.}
}
\value{
none. Output file written into directory of first population comparison.
   (e.g., CEU_YRI directory)
}
\description{
Writes enrichment file for use in Cytoscape EnrichmentMap
}
