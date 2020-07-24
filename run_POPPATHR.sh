#!/bin/bash

# Define population pairs to test
pop_pair_one=CEU_YRI
pop_pair_two=CEU_LWK

# (1) Run gene-set enrichment on defined population comparisons
Rscript get_enrichment.R -p ${pop_pair_one}
Rscript get_enrichment.R -p ${pop_pair_two}

# (2) Run pathway-level coevolution analysis
Rscript get_coevolution.R -p1 ${pop_pair_one} -p2 ${pop_pair_two}

# (3) Run genomic feature integration
Rscript get_properties.R
