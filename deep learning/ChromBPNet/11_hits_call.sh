#!/bin/bash

module load conda
conda activate finemo_gpu

tissue=$1

mkdir -p ATAC_bams/${tissue}/finemo

finemo extract-regions-chrombpnet-h5 -c ATAC_bams/${tissue}/chrombpnet_contribs.counts_scores.h5 -p ATAC_bams/${tissue}/data/peaks_no_blacklist.bed -o ATAC_bams/${tissue}/finemo/${tissue}


finemo call-hits -M pp -r ATAC_bams/${tissue}/finemo/${tissue}.npz --modisco-h5 ATAC_bams/${tissue}/chrombpnet_motifs/${tissue}.modisco_results.h5 -p ATAC_bams/${tissue}/data/peaks_no_blacklist.bed -C BosTaurusARSUCD12V105.chr -o ATAC_bams/${tissue}/finemo


finemo report -r ATAC_bams/${tissue}/finemo/${tissue}.npz -H ATAC_bams/${tissue}/finemo/hits.tsv -p ATAC_bams/${tissue}/data/peaks_no_blacklist.bed -m ATAC_bams/${tissue}/chrombpnet_motifs/${tissue}.modisco_results.h5 -o ATAC_bams/${tissue}/finemo

