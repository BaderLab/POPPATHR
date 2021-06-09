# PopulationPathways #
*The same pathways are present in all humans, but are these genes somehow tuned
differently based on population ancestry?*

This is a project to develop a computational pipeline that determines instances of population-driven SNP-SNP coevolution at the pathway level in an effort to better understand the evolution of human pathways. A wide number of software packages and statistical methods are used and are outlined below (NOTE: the following pipeline is currently carried out between the **CEU** and **YRI** population cohorts)

### Methods summary
| METHOD        | RESULT        |
| ------------- | ------------- |
| PLINK         | ~1.5 million SNPs genotyped in CEU and YRI (HapMap 3)  |
| GSEA          | 56 pathways (19 non-redundant) enriched for population-driven positive selection (i.e., high FST genes)  |
| Linkage disequilibrium | 4 (CEU) / 10 (YRI) within-pathway and 16 (CEU) / 59 (YRI) between-pathway coevolution signals discoered (FDR ≤ 0.2)|

HYPOTHESIS: evolutionary maintenance of pathway-level SNP interactions that influence population fitness

- - - -
## Steps of computational pipeline: ##
Pipeline scripts found [here](https://github.com/rosscm/PopulationPathways/tree/master/bin/R)
```
# within R environment (calls recodeFAM.R, popPCA.R, calcFST.R, SNP2gene.R, setupGSEArun.R, getPathStats.R, LDstatsWPM.R, and LDstatsBPM.R)
> source(runPipeline.R) # runs entire pipeline
```

### 1) PLINK
* PLINK 1.90 for OS X downloaded [here](https://www.cog-genomics.org/plink2)
* Pipeline compatible with **SNP genotyping** data
  * 1,423,541 polymorphic SNPs from the [International HapMap 3 Project](ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/plink_format/)
    * Converted from MAP/PED to BIM/BED/FAM format
  * Genotypes required for two population comparative analyses
    * Original: CEU vs. YRI
    * Replication: CEU vs. LWK
  * Raw PLINK data on Catherine's laptop ~/PopulationPathways/data/HM3_2010-05_phase3
* 3 purposes:
  * Recode PLINK fam file to case/control format by population (`recodeFAM.R`)
  * Perform PCA on tested populations (`popPCA.R`)
  * Calculate SNP-level FST as a marker for population-driven positive selection (`calcFST.R`)

- - - -
### 2) GSEA (Gene-set enrichment analysis)
* Perl scripts from [GenGen](http://gengen.openbioinformatics.org/en/latest/) software tool package

#### Input:
* Pathway annotation file (downloaded from Bader lab [website](http://download.baderlab.org/EM_Genesets/))
  * File used here found in git repo:  [Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt](https://github.com/rosscm/PopulationPathways/blob/master/anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt)
* SNP-level FST statistic file (`calcFST.R`)
  * NOTE: the second column labelled 'CHI2' for purposes of running GenGen GSEA script as it originally expects GWAS *p* values or CHI2, represents FST)
  ```
  > head(markerFST)
  Marker      CHI2
  1  rs1000000  0.100100
  2 rs10000003  0.122500
  3 rs10000007 -0.273025
  4 rs10000010 -0.005300
  5 rs10000011 -0.179690
  6 rs10000015 -0.177080
  ```
* SNP-to-gene mapping file (`SNP2gene.R`)
  * SNPs are mapped to their physically nearest gene
  * GSEA is set to filter out any mapping with SNPs lying >500kb from its nearest gene
  ```
  > head(snp2gene)
  V1           V2 V3
  1 rs12562034 LOC100288069  0
  2  rs7518545 LOC100288069  0
  3  rs9651273         AGRN  0
  4 rs61766341         AGRN  0
  5 rs12134754         AGRN  0
  6 rs11586034         AGRN  0
  ```

#### Output:
* Table containing pathway enrichment results
```
> head(res)
Geneset
1                                                     CELL FATE COMMITMENT%GOBP%GO:0045165
2              SRP-DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE%GOBP%GO:0006614
3                 EUKARYOTIC TRANSLATION TERMINATION%REACTOME DATABASE ID RELEASE 56%72764
4                 EUKARYOTIC TRANSLATION ELONGATION%REACTOME DATABASE ID RELEASE 56%156842
5                                           VIRAL MRNA TRANSLATION%REACTOME%R-HSA-192823.1
6 TRANSMEMBRANE RECEPTOR PROTEIN SERINE/THREONINE KINASE SIGNALING PATHWAY%GOBP%GO:0007178
Size    ES   NES NominalP   FDR   FWER
1  109 0.406 4.667   0.0000 0.009 0.0506
2   83 0.130 4.629   0.0007 0.009 0.0556
3   77 0.114 4.592   0.0025 0.009 0.0623
4   78 0.107 5.024   0.0019 0.010 0.0184
5   75 0.126 4.901   0.0013 0.010 0.0260
6  134 0.372 4.676   0.0000 0.010 0.0492
```

#### Results:
* CEU-YRI: 84 pathways significantly enriched for high-FST genes (FDR ≤ 0.05)
* CEU-LWK: 130 pathways significantly enriched for high-FST genes (FDR ≤ 0.05)
* Replication by both analyses reveals **56 'selection-enriched' pathways** show true evidence of population-driven positive selection
  * In other words, SNPs show significant differentiation in FST between European and African populations at the pathway level
  * *Immune*, *metabolism*, and *development* themes
* Generated binary PLINK files for downstream analysis (`getPathStats.R`)
  * Selection-enriched pathways: 56 (19 non-redundant)
  * Unenriched pathways: 53 (i.e., SNP-level FST values are not significantly differentiated between Europeans and Africans, negative controls)

- - - -
### 3) Linkage disequilibrium (LD)
* Reference [snpStats](http://bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/ld-vignette.pdf) LD vignette
* Inter-chromosomal LD (*r2*) used as a proxy to detect SNP-SNP coevolution at the pathway level
* Two pathway-level interaction models implemented (motivated by recent [BridGE](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006973) method):
  * Within-pathway -- detecting SNP-SNP coevolution within a single pathway
  * Between-pathway -- detecting SNP-SNP coevolution between pairs of pathways (functional pathway clusters)

#### Input:
* SNP genotype and location data per selection enriched and unenriched pathway (`getPathStats.R`)

#### Output:
* LD *r2* value per inter-chromosomal pairwise SNP combination specific to population cohort
  * Within-pathway results
  ```
  > head(dat)
  snp_1 chr_1    pos_1      snp_2 chr_2     pos_2    R.squared   pathway
  1 rs10119465     9 19376144  rs1693547     8 101714394 8.019353e-03 pathway_1
  2 rs10119465     9 19376144  rs1881157     4  99712077 5.271225e-03 pathway_1
  3 rs10119465     9 19376144 rs11684944     2   3633917 1.151087e-02 pathway_1
  4 rs10119465     9 19376144   rs632401     7  73563425 3.049999e-03 pathway_1
  5 rs10119465     9 19376144  rs6545480     2  55460751 2.641132e-02 pathway_1
  6 rs10119465     9 19376144 rs11586570     1  93305371 1.475042e-05 pathway_1
       set pop
  1 Enriched CEU
  2 Enriched CEU
  3 Enriched CEU
  4 Enriched CEU
  5 Enriched CEU
  6 Enriched CEU
  ```
  * Between-pathway results
  ```
  > head(dat)
    snp_1 chr_1     pos_1  gene_1      snp_2 chr_2    pos_2    R.squared
  1   rs166125     5  65201373 ERBB2IP rs10762276    10 71052469 0.0084093756
  2 rs12493155     3 138400332  PIK3CB rs10762276    10 71052469 0.0079737336
  3  rs2131432    17   1297472   YWHAE rs10762276    10 71052469 0.0102617490
  4 rs12076073     1 154944156    SHC1 rs10762276    10 71052469 0.0002102607
  5   rs741334    12 112985328  PTPN11 rs10762276    10 71052469 0.0031810952
  6  rs3756668     5  67596088  PIK3R1 rs10762276    10 71052469 0.0048012803
                                 pathway_pair1
  1 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
  2 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
  3 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
  4 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
  5 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
  6 ALPHA6BETA4INTEGRIN_IOB_ALPHA6BETA4INTEGRIN
                            pathway_pair2       ixn_num      set pop
  1 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  2 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  3 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  4 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  5 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  6 CARBOHYDRATE_TRANSPORT_GOBP_GO_0008643 interaction_1 Enriched CEU
  ```

* Pathway-level FDR statistic representing the significance of coevolution signal (determined via **Kolmogorov-Smirnov** statistic)
  * Within-pathway: LD *r2* distribution per selection-enriched pathway against cumulative LD *r2* distribution of unenriched pathways
  ```
  head(pvals)
  pathway
  1 EMBRYONIC_MORPHOGENESIS
  2 CHROMATIN_SILENCING
  3 NONSENSE_MEDIATED_DECAY_NMD_INDEPENDENT_OF_THE_EXON_JUNCTION_COMPLEX_EJC
  4 NONSENSE_MEDIATED_DECAY_NMD_ENHANCED_BY_THE_EXON_JUNCTION_COMPLEX_EJC
  5 NONSENSE_MEDIATED_DECAY_NMD
  6 CELL_CHEMOTAXIS
   pvals_CEU bonf_CEU fdr_CEU pvals_YRI bonf_YRI fdr_YRI
  1   0.00142   0.0793  0.0793   0.26100    1.000  0.3570
  2   0.00367   0.2058  0.1029   0.80100    1.000  0.8790
  3   0.01135   0.6357  0.1325   0.02480    1.000  0.0555
  4   0.01183   0.6623  0.1325   0.00289    0.162  0.0192
  5   0.01183   0.6623  0.1325   0.00289    0.162  0.0192
  6   0.02210   1.0000  0.1920   0.00253    0.141  0.0192
  ```
    * NOTE: functionally redundant pathways are reported as a single signal of within-pathway coevolution (based on thematic functional grouped as determined by Cytoscape [EnrichmentMap](http://apps.cytoscape.org/apps/enrichmentmap))

  * Between-pathway: LD *r2* distribution per selection-enriched pathway-pathway pair against cumulative LD *r2* distribution of all unenriched pathway-pathway pairs
  ```
  head(pvals)
  interaction.pathway_pair1
  1 ALPHA6BETA4INTEGRIN
  2 EMBRYONIC_MORPHOGENESIS
  3 ALPHA6BETA4INTEGRIN
  4 ALPHA6BETA4INTEGRIN
  5 CELLULAR_RESPONSE_TO_LIGHT
  6 ALPHA6BETA4INTEGRIN
  interaction.pathway_pair2
  1 CELLULAR_RESPONSE_TO_LIGHT
  2 IL4_MEDIATED_SIGNALING_EVENTS
  3 PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS
  4 VIRAL_TRANSLATION_AND_PROTEIN_TARGETING
  5 IMMUNE_CELL_REGULATION                                             
  6 CELL_CHEMOTAXIS
  pvals_CEU bonf_CEU fdr_CEU pvals_YRI bonf_YRI fdr_YRI
  1  0.000379   0.0649  0.0324    0.4730        1   0.691
  2  0.000207   0.0355  0.0324    0.6850        1   0.831
  3  0.002506   0.4285  0.1012    0.2540        1   0.447
  4  0.002959   0.5060  0.1012    0.3830        1   0.583
  5  0.002764   0.4727  0.1012    0.8360        1   0.929
  6  0.013680   1.0000  0.1533    0.0268        1   0.112
  ```

#### Results:
* Within-pathway: 4 (CEU) and 10 (YRI) non-redundant significant coevolution signals discovered
* Between-pathway: 16 (CEU) and 59 (YRI) significant coevolution signals discovered
* Primarily **immune-associated** themes

**SUMMARY**
* Proof-of-concept pipeline successfully discovers instances of coevolution among pathways influencing population fitness
  * Pathways demonstrate population-driven signals of positive selection (FST) between Europeans and Africans
  * LD maintains inter-chromosomal SNP-SNP interactions within and between pathways enriched for positive selection
