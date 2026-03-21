#!/bin/bash

mark=$1
tissue=$2
module load R
Rscript prepare_neg_pos_sets.R ${mark} ${tissue} negative
