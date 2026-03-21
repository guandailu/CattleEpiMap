#!/bin/bash


tissue=$1
#chr=$2

module load conda 
conda activate chrombpnet
module load apptainer
module load cuda/11.7.1 
export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"
mkdir -p ATAC_bams/${tissue}/SNP_scores ATAC_bams/${tissue}/Cattle_SNPs/
#cp -r Cattle_SNPs/chr${chr}.list ATAC_bams/${tissue}/Cattle_SNPs/chr${chr}.list

#head -n -30 Cattle_SNPs/chr${chr}.list > ATAC_bams/${tissue}/Cattle_SNPs/chr${chr}.list

#singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet snp_score -snps /mnt/Cattle_SNPs/chr${chr}.list -m /mnt/chrombpnet_model/models/chrombpnet.h5 -g /mnt/data/ref.fa -op /mnt/SNP_scores/chr${chr}

chr_subfile=$2
output_file=$(echo ${chr_subfile} | sed 's/.*\///g')


chrombpnet snp_score -snps ${chr_subfile} -m ATAC_bams/${tissue}/chrombpnet_model/models/chrombpnet_nobias.h5 -g ATAC_bams/${tissue}/data/ref.fa -op ATAC_bams/${tissue}/SNP_scores/${output_file}
