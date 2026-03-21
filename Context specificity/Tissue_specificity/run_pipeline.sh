#!/bin/bash -l

set -e

state=$1
path="./"
mkdir -p $state
module load bedtools2
cat ${path}/03_chromatin_ann/FourMarks/ChromatinStates/*_States.bed  | awk -v e=$state '{if ($4 == e) print}' | bedtools sort -i - | bedtools merge -i - -d 200 > ${state}/consensus_state.bed

# create consensus element regions
mkdir -p ${state}/consensus_state
for i in `ls -1 ${path}/03_chromatin_ann/FourMarks/ChromatinStates/*_States.bed | sed 's/_/\t/g' | cut -f3 | sed 's/.*\///g;s/AW//g;s/BW//g;s/btr//g;s/ctl//g' | sort | uniq`; 
do
    cat ${path}/03_chromatin_ann/FourMarks/ChromatinStates/${i}*_States.bed | awk -v e=$state '{if ($4 == e) print}' | bedtools sort -i - | bedtools merge -i - -d 200 | awk '{print $0"\t"1}' | bedtools map -a ${state}/consensus_state.bed -b - -c 4 -null 0 -F 0.5 > ${state}/consensus_state/${i}.bed
done

## manually prepare tissues_domain.txt, the first column is tissue name and the second is domain name
#identify tissue specific element
for t in `ls -1 ${path}/03_chromatin_ann/FourMarks/ChromatinStates/*_States.bed | sed 's/_/\t/g' | cut -f3 | sed 's/.*\///g;s/AW//g;s/BW//g;s/btr//g;s/ctl//g' | sort | uniq`;
do
    num=$(cat tissues_domain.txt | awk -v t=$t '{if ($1 == t) print $2}' | wc -l)
    if [[ $num -eq 1 ]]
    then
        rest_tiss=$(cat tissues_domain.txt | awk -v t=$t -v e=$state '{if ($1 != t) print e"/consensus_state/"$1"*.bed"}' | tr '\n' ' ')
    else
        domain=$(cat tissues_domain.txt | awk -v t=$t '{if ($1 == t) print $2}' | uniq)
        rest_tiss=$(cat tissues_domain.txt | awk -v d=$domain -v e=$state '{if ($2 == d) print e"/consensus_state/"$1"*.bed"}' | tr '\n' ' ')
    fi
    cat ${state}/consensus_state/${t}.bed | awk '{if ($4 > 0) print}' | bedtools sort -i - | bedtools intersect -a - -b <(cat ${rest_tiss} | awk '{if ($4 >0 ) print}' | bedtools sort -i - ) -v > ${state}/consensus_state/${t}_specific_${state}.txt
done


module load R
Rscript tissue_specific_RE_analysis.R $state
module unload R

source ~/miniconda3/bin/activate homer
mkdir -p Temp ${state}/Tisssue_specific_HOMER
for i in ${state}/consensus_state/*_specific_${state}.txt
do
    n=$(echo $i | sed "s/.*consensus_state\///g;s/_specific_.*.txt//g")
    tmpbed=$(mktemp Temp/bed.XXXXXXXXXXXX)
    cat $i | awk '{if ($4 == 1) print}' > ${tmpbed}
    findMotifsGenome.pl ${tmpbed} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.fa ${state}/Tisssue_specific_HOMER/${n} -size 200 -len 8 -p 12
    rm  -rf ${tmpbed}
done

