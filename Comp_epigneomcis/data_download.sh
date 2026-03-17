#!/bin/bash -l

set -e

wget -m -nd -A "impute_*_H3K27ac.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
wget -m -nd -A "impute_*_H3K27me3.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
wget -m -nd -A "impute_*_H3K4me1.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
wget -m -nd -A "impute_*_H3K4me3.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
wget -m -nd -A "impute_*_ATAC*.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
wget -m -nd -A "impute_*_CTCF*.bigWig" https://epigenome.wustl.edu/epimap/data/imputed/
