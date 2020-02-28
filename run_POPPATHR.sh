#!/bin/bash

# Run gene-set enrichment on two population comparisons
pop_pair_one=CEU_YRI
Rscript run_enrichment.R -p ${pop_pair_one} --SET_PERM 100

pop_pair_two=CEU_LWK
Rscript get_enrichment.R -p ${pop_pair_two} --SET_PERM 100

# Run coevolution analysis
Rscript get_coevolution.R -p1 CEU_YRI -p2 CEU_LWK

# Run gene property integration
Rscript get_properties.R
