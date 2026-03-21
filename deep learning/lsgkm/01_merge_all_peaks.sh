#!/bin/bash

mark=$1
mkdir -p data/${mark}

module load bedtools2

cat 01_data_preprocess/07_peak_called/${mark}_*_Peaks.bed  | awk -v OFS="\t" '{if ($1 >= 1 && $1 <= 29) print "chr"$1, $2, $3}' | bedtools sort -i - | bedtools merge -i - > data/all_samples_${mark}_peaks.bed

