#!/bin/bash


tissue=$1
#chr=$2

module load conda 
conda activate chrombpnet
module load apptainer
module load cuda/11.7.1 
module avail cudnn
export APPTAINER_CACHEDIR="23_deep_learning/chrombpnet/ATAC_bams/${tissue}/APPTAINER_CACHEDIR"
mkdir -p ATAC_bams/${tissue}/ASE_scores_bias_model
#cp -r Cattle_SNPs/chr${chr}.list ATAC_bams/${tissue}/Cattle_SNPs/chr${chr}.list

#head -n -30 Cattle_SNPs/chr${chr}.list > ATAC_bams/${tissue}/Cattle_SNPs/chr${chr}.list

#singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet snp_score -snps /mnt/Cattle_SNPs/chr${chr}.list -m /mnt/chrombpnet_model/models/chrombpnet.h5 -g /mnt/data/ref.fa -op /mnt/SNP_scores/chr${chr}

chr_subfile=$2
output_file=$(echo ${chr_subfile} | sed 's/.*\///g')


#chrombpnet snp_score -snps Cattle_SNPs/chr${chr}_split.0000000028 -m ATAC_bams/${tissue}/chrombpnet_model/models/chrombpnet.h5 -g ATAC_bams/${tissue}/data/ref.fa -op ATAC_bams/${tissue}/SNP_scores/chr${chr}_split.0000000028


#chrombpnet snp_score -snps ${chr_subfile} -m ATAC_bams/${tissue}/chrombpnet_model/models/chrombpnet_nobias.h5 -g ATAC_bams/${tissue}/data/ref.fa -op ATAC_bams/${tissue}/ASE_scores/${output_file}
cp $chr_subfile ATAC_bams/${tissue}/.
singularity  exec --nv -e --no-mount hostfs --bind ATAC_bams/${tissue}/:/mnt docker://kundajelab/chrombpnet:latest chrombpnet snp_score -snps /mnt/${chr_subfile} -m /mnt/chrombpnet_model/models/chrombpnet.h5 -g /mnt/data/ref.fa -op /mnt/ASE_scores_bias_model/${output_file}

rm -rf ATAC_bams/${tissue}/${chr_subfile}
