#!/bin/bash

tissue=$1

module load conda 
conda activate chrombpnet
module load apptainer

export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

mkdir -p ATAC_bams/${tissue}/pred_bw
singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet pred_bw -bm ATAC_bams/${tissue}/chrombpnet_model/models/bias_model_scaled.h5 -cm ATAC_bams/${tissue}/chrombpnet_model/models/chrombpnet.h5 -cmb ATAC_bams/${tissue}/chrombpnet_model/models/chrombpnet_nobias.h5 -r ATAC_bams/${tissue}/data/peaks_no_blacklist.bed -g ATAC_bams/${tissue}/data/ref.fa -c ATAC_bams/${tissue}/data/ref.chrom.sizes -op ATAC_bams/${tissue}/pred_bw/${tissue}_peaks 

