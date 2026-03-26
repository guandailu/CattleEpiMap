#!/bin/bash

set -e

echo -e "************ Parse parameters ************ "
sample=$1
echo -e "Processing ${sample}...\n\n"
mkdir -p Temp

function read_type(){
    data_path1="ChIPQC/Raw_Reads"
    data_path2="PublicChIP/Raw_Reads"
    r2=${sample}_R2.fq.gz
    if [[ -f ${data_path1}/${r2} ]] && [[ ! -f ${data_path2}/${r2} ]]
    then
        echo -e "Raw data is in the folder: ${data_path1}\n"
        fastq_r1=${data_path1}/${sample}_R1.fq.gz
        fastq_r2=${data_path1}/${sample}_R2.fq.gz
        seq_type="paired"
    elif [[ ! -f ${data_path1}/${r2} ]] && [[ -f ${data_path2}/${r2} ]]
    then
       echo -e "Raw data is in the folder: ${data_path2}\n"
       fastq_r1=${data_path2}/${sample}_R1.fq.gz
       fastq_r2=${data_path2}/${sample}_R2.fq.gz
       seq_type="paired"
    elif [[ ! -f ${data_path1}/${r2} ]] && [[ ! -f ${data_path2}/${r2} ]]
    then
        echo -e "No second end found, fastq file should be single-end...\n"
        rr=${sample}.fq.gz
        if [[ -f ${data_path1}/${rr} ]] && [[ ! -f ${data_path2}/${rr} ]]
        then
            fastq_r1=${data_path1}/${rr}
            seq_type="single"
        elif [[ ! -f ${data_path1}/${rr} ]] && [[ -f ${data_path2}/${rr} ]]
        then
            fastq_r1=${data_path2}/${rr}
            seq_type="single"
        else
            echo -e "Error: No read files for ${sample}, please double-check ${sample} raw data...\n"
            exit 1
        fi
    fi
    echo -e "${sample} is ${seq_type}-end...\n\n\n"
}

function run_fastqc(){
    mkdir -p 00_fastqc_raw/${sample}
    if [ ${seq_type} == "paired" ]
    then
        fastqc -t 8 -o 00_fastqc_raw/${sample} ${fastq_r1} ${fastq_r2}
    elif [ ${seq_type} == "single" ]
    then
        fastqc -t 8 -o 00_fastqc_raw/${sample} ${fastq_r1}
    else 
        echo -e "Error: fastqc raw reads for ${sample}...\n"
        exit 1
    fi
}

function trimming_reads(){
    mkdir -p 01_trimmed_reads/${sample}
    if [ ${seq_type} == "paired" ]
    then
        #trim_galore --paired --trim-n --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --length 30 -q 30 --gzip --trim-n -j 8 -o 01_trimmed_reads --basename ${sample} --fastqc --fastqc_args "--outdir 01_trimmed_reads/${sample} -t 8"  ${fastq_r1} ${fastq_r2}
        trim_galore --paired --trim-n --length 30 -q 30 --gzip --trim-n -j 8 -o 01_trimmed_reads --basename ${sample} --fastqc --fastqc_args "--outdir 01_trimmed_reads/${sample} -t 8"  ${fastq_r1} ${fastq_r2}
  unzip -o -d 01_trimmed_reads/${sample} 01_trimmed_reads/${sample}/${sample}_val_1_fastqc.zip
        num_r1=$(cat 01_trimmed_reads/${sample}/${sample}_val_1_fastqc/fastqc_data.txt | awk '{if (NR==7) print $3}')
        echo -e "The number of reads aftering trimming is: ${num_r1}\n"
    elif [ ${seq_type} == "single" ]
    then
        #trim_galore -q 30 --trim-n --clip_R1 5 --three_prime_clip_R1 5 --length 30 --gzip --trim-n -j 8 -o 01_trimmed_reads --basename ${sample} --fastqc --fastqc_args "--outdir 01_trimmed_reads/${sample} -t 8" ${fastq_r1}
  trim_galore -q 30 --trim-n --length 30 --gzip --trim-n -j 8 -o 01_trimmed_reads --basename ${sample} --fastqc --fastqc_args "--outdir 01_trimmed_reads/${sample} -t 8" ${fastq_r1}
        unzip -o -d 01_trimmed_reads/${sample} 01_trimmed_reads/${sample}/${sample}_trimmed_fastqc.zip
        num_r1=$(cat 01_trimmed_reads/${sample}/${sample}_trimmed_fastqc/fastqc_data.txt | awk '{if (NR==7) print $3}')
        echo -e "The number of reads aftering trimming is: ${num_r1}\n"
    else
        echo -e "Error: Sequencing type {single|paired} is not correct ...\n"
        exit 1
    fi
}

function mapping_reads(){
    rg="@RG\\tID:${sample}\\tSM:${sample}"
    fasta="/group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.fa"
    if [[ ${sample} =~ "ATAC" ]] || [[ ${seq_type} == "paired" ]]
    then
        bwa mem -M -t 8 -R ${rg} ${fasta} ${input_fq1} ${input_fq2} > Temp/${sample}.aligned.sam
        samtools sort -O BAM -@ 8 Temp/${sample}.aligned.sam > Temp/${sample}.sorted.bam
    else
        bwa mem -M -t 8 -R ${rg} ${fasta} ${input_fq} > Temp/${sample}.aligned.sam
        samtools sort -O BAM -@ 8  Temp/${sample}.aligned.sam > Temp/${sample}.sorted.bam
    fi
    samtools index Temp/${sample}.sorted.bam
}

function filter_alignment(){
    mkdir -p 03_filtered_alignments
    samtools idxstats Temp/${sample}.sorted.bam | cut -f 1 | grep -E "^[0-9]" | xargs samtools view -b Temp/${sample}.sorted.bam > Temp/${sample}.intermediate.bam      
    samtools view -h -F 1804 -q 30 -@ 8 Temp/${sample}.intermediate.bam | grep -v XA:Z | grep -v SA:Z | samtools sort -@ 8 -n - > Temp/${sample}.filtered.bam
    samtools fixmate -@ 8 -m Temp/${sample}.filtered.bam Temp/${sample}.fixmate.bam
    samtools sort -@ 8 -o Temp/${sample}.sorted2.bam Temp/${sample}.fixmate.bam
    samtools markdup -@ 8 -O BAM -r --write-index Temp/${sample}.sorted2.bam 03_filtered_alignments/${sample}.dedup.bam
}

function spp_stats(){
    mkdir -p 06_metrics
    Rscript /home/dguan/bin/FunctionalAnnotation/Scripts/run_spp.R -c=05_mappability_filter/${sample}.bam -rf -out=06_metrics/${sample}.spp_stats.txt -p=8 -s=0:2:400 -savp=06_metrics/${sample}.Cross_Correlation.pdf -tmpdir=Temp
    samtools view -c -F 260 05_mappability_filter/${sample}.bam > 06_metrics/${sample}.num_reads.txt
}


echo -e "************ Determing reads type ************ \n\n\n"
read_type


echo -e "************ Trimming reads ************ \n"
module load trimgalore/0.6.6
trimming_reads
module unload trimgalore/0.6.6
if [ ${seq_type} == "paired" ]
then
    input_fq1="../01_data_preprocess_2/01_trimmed_reads/${sample}_val_1.fq.gz"
    input_fq2="../01_data_preprocess_2/01_trimmed_reads/${sample}_val_2.fq.gz"
elif [ ${seq_type} == "single" ]
then
    input_fq="../01_data_preprocess_2/01_trimmed_reads/${sample}_trimmed.fq.gz"
fi

echo -e "************ Mapping reads against the reference ************ \n\n\n"
module load bwa 
module load samtools 
mapping_reads
module unload bwa


echo -e "************ Filtering alignment ************ \n\n\n"
module load bedtools2/2.30.0
filter_alignment
module unload bedtools2/2.30.0


echo -e "************ Filtering by data mappability ************ \n\n\n"
mkdir -p 05_mappability_filter
module load bedtools2
MAP36="BosTaurusARSUCD12V105_mappability_ReadLength36.k36.bed"
bedtools intersect -abam 03_filtered_alignments/${sample}.dedup.bam -b ${MAP36} -wa > 05_mappability_filter/${sample}.bam
bedtools bamtobed -i 05_mappability_filter/${sample}.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > 05_mappability_filter/${sample}.tagAlign.gz
module unload bedtools2

echo -e "************ Collecting data metrics ************ \n\n\n"
module load R/4.2.3
module load fastqc/0.11.9
if [[ ! ${sample} =~ Input* ]]
then
    spp_stats
fi
module unload R/4.2.3
module unload fastqc/0.11.9
module unload samtools
