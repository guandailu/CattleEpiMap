#!/bin/bash

tissue=$1

#### Training bias model
module load conda 
conda activate chrombpnet
module load apptainer
module load cuda/11.7.1 


export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"


singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet contribs_bw -m /mnt/chrombpnet_model/models/chrombpnet_nobias.h5 -r /mnt/data/peaks_no_blacklist.bed -g /mnt/data/ref.fa -c /mnt/data/ref.chrom.sizes -op /mnt/chrombpnet_contribs



