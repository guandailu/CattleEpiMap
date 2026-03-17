#!/bin/bash

mark=$1
tissue=$2

#cd ${tissue}
mkdir -p data/${mark}/${tissue}

module load R
Rscript prepare_neg_pos_sets.R ${mark} ${tissue} positive

gkmtrain -T 16 data/${mark}/${tissue}/${mark}_${tissue}_Peaks.fa data/all_samples_${mark}_peaks.neg1x.fa data/${mark}/${tissue}/${mark}_${tissue}


