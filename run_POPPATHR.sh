#!/bin/bash

# Define population pairs to test
pop_pair_one=CEU_YRI
pop_pair_two=CEU_LWK

# Run gene-set enrichment on defined populations
Rscript get_enrichment.R -p ${pop_pair_one} --SET_PERM 100
Rscript get_enrichment.R -p ${pop_pair_two} --SET_PERM 100

# Run coevolution analysis
Rscript get_coevolution.R -p1 ${pop_pair_one} -p2 ${pop_pair_two}

# Run gene property integration
Rscript get_properties.R
