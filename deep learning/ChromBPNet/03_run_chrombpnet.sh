#!/bin/bash -l

tissue=$1
rm -rf ATAC_bams/${tissue}/bias_model
rm -rf ATAC_bams/${tissue}/chrombpnet_model/
#### Training bias model
module load conda 
conda activate chrombpnet
module load apptainer
module load cuda/11.7.1 
mkdir -p ATAC_bams/${tissue}/APPTAINER_CACHEDIR
export APPTAINER_CACHEDIR="/group/zhougrp2/dguan/BovineFAANG/23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"
#### Bias model training
singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet bias pipeline \
        -ibam /mnt/data/merged.bam \
        -d "ATAC" \
        -g /mnt/data/ref.fa \
        -c /mnt/data/ref.chrom.sizes \
        -p /mnt/data/peaks_no_blacklist.bed \
        -n /mnt/data/output_negatives.bed \
        -fl /mnt/data/splits/fold_0.json \
        -b 0.5 \
        -o /mnt/data/bias_model/ \
        -fp ${tissue} \


#### Training bias-factorized ChromBPNet
singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet pipeline -ibam /mnt/data/merged.bam -d ATAC -g /mnt/data/ref.fa -c /mnt/data/ref.chrom.sizes -p /mnt/data/peaks_no_blacklist.bed -n /mnt/data/output_negatives.bed -fl /mnt/data/splits/fold_0.json -b /mnt/data/bias_model/models/${tissue}_bias.h5 -o /mnt/chrombpnet_model/
