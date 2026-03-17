#!/bin/bash

module load conda 
conda activate chrombpnet
#module load apptainer
#export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"
tissue=$1
mkdir -p ATAC_bams/${tissue}/chrombpnet_motifs
modisco motifs -i ATAC_bams/${tissue}/chrombpnet_contribs.counts_scores.h5 -n 1000000 -o ATAC_bams/${tissue}/chrombpnet_motifs/${tissue}.modisco_results.h5

modisco report -i ATAC_bams/${tissue}/chrombpnet_motifs/${tissue}.modisco_results.h5 -o ATAC_bams/${tissue}/chrombpnet_motifs/ -m /group/zhougrp4/dguan/BovineFAANG/23_deep_learning/chrombpnet/JASPAR2020_CORE_non-redundant_pfms_meme.txt

