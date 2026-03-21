#!/bin/bash -l

set -e

path="./"

module load bedtools2
mkdir -p tissue_specificity
se_files=$(cat ROSE/tissues.list | grep -v "Liver" | awk '{print "output/"$1"/enhancer_SuperStitched.table.txt"}' | tr '\n' ' ')
cat $se_files | grep -v "^#" | grep -v "REGION_ID" | awk -v OFS="\t" '{print $2,$3,$4}' | sed 's/chr//g' | bedtools sort -i - | bedtools merge -i - -d 10000 > tissue_specificity/consensus_SE.bed

mkdir -p tissue_specificity/consensus_SE
for i in `cat ROSE/tissues.list | grep -v "Liver"`; 
do
    cat ${path}/04_super_enhacner/output/${i}/enhancer_SuperStitched.table.txt | grep -v "^#" | grep -v "REGION_ID" | cut -f2-4 | sed 's/chr//g' | bedtools sort -i - | awk '{print $0"\t"1}' | bedtools map -a tissue_specificity/consensus_SE.bed -b - -c 4 -null 0 -F 0.5 | awk '{if ($4 > 0) print}' > tissue_specificity/consensus_SE/${i}_SEs.bed
done

## manually prepare tissues_domain.txt, the first column is tissue name and the second is domain name
for t in `cat ROSE/tissues.list | grep -v "Liver"`;
do
    num=$(cat tissues_domain.txt | grep -v "Liver" | awk -v t=$t '{if ($1 == t) print $2}' | wc -l)
    if [[ $num -eq 1 ]]
    then
        rest_tiss=$(cat tissues_domain.txt | grep -v "Liver" | awk -v t=$t '{if ($1 != t) print "tissue_specificity/consensus_SE/"$1"_SEs.bed"}' | tr '\n' ' ')
    else
        domain=$(cat tissues_domain.txt | grep -v "Liver" | awk -v t=$t '{if ($1 == t) print $2}' | uniq)
        rest_tiss=$(cat tissues_domain.txt | grep -v "Liver" | awk -v d=$domain '{if ($2 == d) print "tissue_specificity/consensus_SE/"$1"_SEs.bed"}' | tr '\n' ' ')
    fi
    cat tissue_specificity/consensus_SE//${t}_SEs.bed | awk '{if ($4 == 1) print}' | bedtools sort -i - | bedtools intersect -a - -b <(cat ${rest_tiss} | awk '{if ($4 == 1) print}' | bedtools sort -i - ) -v > tissue_specificity/consensus_SE/${t}_specific_SEs.txt
done

source ~/miniconda3/bin/activate homer
mkdir -p Temp tissue_specificity/Tisssue_specific_HOMER
for i in tissue_specificity/consensus_SE/*_specific_SEs.txt
do
    n=$(echo $i | sed "s/.*consensus_SE\///g;s/_specific_SEs.txt//g")
    tmpbed=$(mktemp Temp/bed.XXXXXXXXXXXX)
    cat $i | awk '{if ($4 == 1) print}' > ${tmpbed}
    findMotifsGenome.pl ${tmpbed} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.fa tissue_specificity/Tisssue_specific_HOMER/${n} -size 200 -len 8 -p 12
    rm  -rf ${tmpbed}
done

#
