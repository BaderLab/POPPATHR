# ~/PopulationPathways/anno/baderlab

source("~/PopulationPathways/bin/R/post_plots/readPathways.R")
library(cogena) #write gmt list to file

setwd("~/PopulationPathways/anno/baderlab")
annoF <- "Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol"
paths <- readPathways(sprintf("%s.gmt", annoF), MIN_SIZE=0, MAX_SIZE=100000)

paths[["CELL CHEMOTAXIS"]] <- unique(c(paths[["CELL CHEMOTAXIS"]],
                                       paths[["CHEMOTAXIS"]],
                                       paths[["TAXIS"]]))
paths[["CELLULAR RESPONSE TO LIGHT"]] <- unique(c(paths[["CELLULAR RESPONSE TO LIGHT STIMULUS"]],
                                                  paths[["CELLULAR RESPONSE TO RADIATION"]]))
paths[["GROWTH FACTOR AND BMP STIMULUS"]] <- unique(c(paths[["BMP SIGNALING PATHWAY"]],
                                                      paths[["CELLULAR RESPONSE TO BMP STIMULUS"]],
                                                      paths[["CELLULAR RESPONSE TO GROWTH FACTOR STIMULUS"]],
                                                      paths[["REGULATION OF PATHWAY-RESTRICTED SMAD PROTEIN PHOSPHORYLATION"]],
                                                      paths[["RESPONSE TO BMP"]],
                                                      paths[["RESPONSE TO GROWTH FACTOR"]],
                                                      paths[["SIGNALING BY BMP"]],
                                                      paths[["SMAD PROTEIN SIGNAL TRANSDUCTION"]],
                                                      paths[["TRANSMEMBRANE RECEPTOR PROTEIN SERINE/THREONINE KINASE SIGNALING PATHWAY"]]))
paths[["IMMUNE CELL REGULATION"]] <- unique(c(paths[["REGULATION OF HEMOPOIESIS"]],
                                              paths[["REGULATION OF MYELOID CELL DIFFERENTIATION"]],
                                              paths[["REGULATION OF MYELOID LEUKOCYTE DIFFERENTIATION"]]))
paths[["REGULATION OF CELL MIGRATION"]] <- unique(c(paths[["POSITIVE REGULATION OF CELL MIGRATION"]],
                                                    paths[["POSITIVE REGULATION OF CELL MOTILITY"]]))
paths[["REGULATION OF LIPID METABOLISM"]] <- unique(c(paths[["POSITIVE REGULATION OF LEUKOCYTE MIGRATION"]],
                                                      paths[["POSITIVE REGULATION OF LIPID METABOLIC PROCESS"]],
                                                      paths[["REGULATION OF LIPID METABOLIC PROCESS"]]))
paths[["VIRAL TRANSLATION AND PROTEIN TARGETING"]] <- unique(c(paths[["3 -UTR-MEDIATED TRANSLATIONAL REGULATION"]],
                                                               paths[["CAP-DEPENDENT TRANSLATION INITIATION"]],
                                                               paths[["COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE"]],
                                                               paths[["ESTABLISHMENT OF PROTEIN LOCALIZATION TO ORGANELLE"]],
                                                               paths[["EUKARYOTIC TRANSLATION ELONGATION"]],
                                                               paths[["EUKARYOTIC TRANSLATION INITIATION"]],
                                                               paths[["EUKARYOTIC TRANSLATION TERMINATION"]],
                                                               paths[["FORMATION OF A POOL OF FREE 40S SUBUNITS"]],
                                                               paths[["GTP HYDROLYSIS AND JOINING OF THE 60S RIBOSOMAL SUBUNIT"]],
                                                               paths[["L13A-MEDIATED TRANSLATIONAL SILENCING OF CERULOPLASMIN EXPRESSION"]],
                                                               paths[["NONSENSE MEDIATED DECAY (NMD) ENHANCED BY THE EXON JUNCTION COMPLEX (EJC)"]],
                                                               paths[["NONSENSE MEDIATED DECAY (NMD) INDEPENDENT OF THE EXON JUNCTION COMPLEX (EJC)"]],
                                                               paths[["NONSENSE-MEDIATED DECAY (NMD)"]],
                                                               paths[["PEPTIDE CHAIN ELONGATION"]],
                                                               paths[["PROTEIN IMPORT INTO NUCLEUS"]],
                                                               paths[["PROTEIN TARGETING TO ER"]],
                                                               paths[["PROTEIN TARGETING TO NUCLEUS"]],
                                                               paths[["SELENOCYSTEINE SYNTHESIS"]],
                                                               paths[["SINGLE-ORGANISM NUCLEAR IMPORT"]],
                                                               paths[["SRP-DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE"]],
                                                               paths[["VIRAL MRNA TRANSLATION"]]))
paths[["WNT CALCIUM SIGNALLING"]] <- unique(c(paths[["CA2+ PATHWAY"]],
                                              paths[["WNT SIGNALING PATHWAY"]]))
cat(sprintf("Writing out new gmt file %s_EM_clusters.gmt...", annoF))
gmtlist2file(paths, sprintf("%s_EM_clusters.gmt", annoF))
cat("done.\n")
