#!/bin/bash
# Script to download genome and pathway annotation files used in manuscript
SCRIPT_DIR=$(pwd)
ANNO_DIR=$SCRIPT_DIR/anno
cd $ANNO_DIR

date="April_24_2016"
organism="Human"
key_type="symbol"
path_home="http://download.baderlab.org/EM_Genesets/"
path_file="${organism}_GOBP_AllPathways_no_GO_iea_${date}_${key_type}.gmt"

wget ${path_home}/${date}/${organism}/${key_type}/${path_file}
