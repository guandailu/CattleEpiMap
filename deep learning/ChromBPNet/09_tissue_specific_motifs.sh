#!/bin/bash

tissue=$1

#### Training bias model
module load conda 
conda activate chrombpnet
module load apptainer
module load cuda/11.7.1 


export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"


mkdir -p ATAC_bams/${tissue}/tissue_specific_motifs

singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet contribs_bw -m /mnt/chrombpnet_model/models/chrombpnet_nobias.h5 -r /mnt/data/tissue_specific_peaks.bed -g /mnt/data/ref.fa -c /mnt/data/ref.chrom.sizes -op /mnt/tissue_specific_motifs/${tissue}.chrombpnet_contribs




modisco motifs -i ATAC_bams/${tissue}/tissue_specific_motifs/${tissue}.chrombpnet_contribs.counts_scores.h5 -n 1000000 -o ATAC_bams/${tissue}/tissue_specific_motifs/${tissue}.modisco_results.h5

modisco report -i ATAC_bams/${tissue}/tissue_specific_motifs/${tissue}.modisco_results.h5 -o ATAC_bams/${tissue}/tissue_specific_motifs/chrombpnet_motifs/ -m /group/zhougrp4/dguan/BovineFAANG/23_deep_learning/chrombpnet/JASPAR2020_CORE_non-redundant_pfms_meme.txt

