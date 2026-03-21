#!/bin/bash

module load py-deeptools/3.5.3 

tissue=$1
bamCoverage --bam ATAC_bams/${tissue}/data/merged.withchr.bam -o ATAC_bams/${tissue}/data/merged.withchr.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2489385779 --ignoreForNormalization chrX --extendReads -p 12

cat ATAC_bams/${tissue}/data/peaks_no_blacklist.bed | awk -v OFS="\t" '{print $1,$2,$3,$1"_"$2"_"$3}' > ATAC_bams/${tissue}/data/peaks_no_blacklist.observed.bed

bigWigAverageOverBed ATAC_bams/${tissue}/data/merged.withchr.bw ATAC_bams/${tissue}/data/peaks_no_blacklist.observed.bed ATAC_bams/${tissue}/data/peaks_no_blacklist.observed.out

bigWigAverageOverBed ATAC_bams/${tissue}/pred_bw/Lung_peaks_chrombpnet_nobias.bw ATAC_bams/${tissue}/data/peaks_no_blacklist.observed.bed ATAC_bams/${tissue}/pred_bw/peaks_no_blacklist.predicted.bed 

