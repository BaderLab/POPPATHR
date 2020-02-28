
## STEP 3: Getting gene properties for selection-enriched pathways ##
message("\n**Cross-referencing selection-enriched interactions with BioGRID.\n")
bgridDir <- sprintf("%s/biogrid", outDir)
if (!file.exists(bgridDir)) dir.create(bgridDir)

enrich_file <- sprintf("%s/genes_enrich_unique.txt", enrich_folder)
getBiogrid(bgridF=bgridF, refgene_file=enrich_file, outDir=bgridDir)


message("\n**Grabbing selection-enriched gene properties.\n")
propDir <- sprintf("%s/geneProp", outDir)
if (!file.exists(propDir)) dir.create(propDir)

geneDir <- enrich_folder
refgene_file  <- sprintf("%s/genes_leadingEdge_unique_hc.txt", geneDir)
clustF <- sprintf("%s/gseaStat_per_hc_pathway_updated.txt", geneDir)
outDir=propDir
