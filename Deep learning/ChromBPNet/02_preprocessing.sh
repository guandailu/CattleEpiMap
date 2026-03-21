#!/bin/bash -l

tissue=$1
path="./"

rep1=$(cat ${path}/ATAC_samples.list | grep -E "_${tissue}_" | head -n 1)
rep2=$(cat ${path}/ATAC_samples.list | grep -E "_${tissue}_" | tail -n 1)

#### merge bam files by tissue
bam1="./${rep1}.bam"
bam2="./${rep2}.bam"
rm -rf ATAC_bams/${tissue}
mkdir -p ATAC_bams/${tissue}/data
module load samtools
samtools merge -@12 -f ATAC_bams/${tissue}/merged_unsorted.bam ${bam1}  ${bam2}
samtools sort -@12 ATAC_bams/${tissue}/merged_unsorted.bam -o ATAC_bams/${tissue}/data/merged.nochr.bam
samtools view -H  ATAC_bams/${tissue}/data/merged.nochr.bam | awk '/^@SQ/ {sub(/SN:/, "SN:chr");} {print}' | samtools reheader - ATAC_bams/${tissue}/data/merged.nochr.bam > ATAC_bams/${tissue}/data/merged.withchr.bam
samtools index -@12 ATAC_bams/${tissue}/data/merged.withchr.bam
samtools view -b ATAC_bams/${tissue}/data/merged.withchr.bam chr{1..29} > ATAC_bams/${tissue}/data/merged.bam
samtools index -@12 ATAC_bams/${tissue}/data/merged.bam
sleep 3
rm -rf ATAC_bams/${tissue}/merged_unsorted.bam
rm -rf ATAC_bams/${tissue}/merged.nochr.bam
#### run idr to get consensus peaks for two replicates
mkdir -p ATAC_bams/${tissue}/idr_output/
mark="ATAC"
mark_type="narrowPeak"
module load conda
conda activate idr
idr --samples ${path}/ATAC_narrow_peaks/${rep1}_Peaks.bed.gz ${path}/ATAC_narrow_peaks/${rep2}_Peaks.bed.gz \
    --input-file-type ${mark_type} \
    --output-file ATAC_bams/${tissue}/idr_output/${mark}_${tissue}.idr.txt \
    --plot \
    --log-output-file ATAC_bams/${tissue}/idr_output/${mark}_${tissue}.log
module load bedtools2
bedtools intersect -wo -a ATAC_bams/${tissue}/idr_output/${mark}_${tissue}.idr.txt -b <(zcat ${path}/ATAC_narrow_peaks/${rep1}_Peaks.bed.gz) | awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$23-$22; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}' | cut -f 1-20 | sort -u | bedtools intersect -wo -a stdin -b <(zcat ${path}/ATAC_narrow_peaks/${rep2}_Peaks.bed.gz) |  awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$23-$22; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}' | cut -f 1-20 | sort -u | awk '{if ($1 >=1 && $1 <= 29) print "chr"$0}' > ATAC_bams/${tissue}/idr_output/${mark}_${tissue}_overlap_peak.bed
module unload conda

#### remove black list from peaks
cat refGenome/BosTaurusARSUCD12V105.mappability_ReadLength36.blacklist.bed | awk '{print "chr"$0}' | bedtools sort -i - | bedtools merge -i -  > ATAC_bams/${tissue}/data/blacklist.bed
blacklist="ATAC_bams/${tissue}/data/blacklist.bed"
cat refGenome/BosTaurusARSUCD12V105.chr.fa.genome | awk '{if ($1 >=1 && $1 <= 29) print "chr"$0}' > ATAC_bams/${tissue}/data/ref.chrom.sizes
chrsize="ATAC_bams/${tissue}/data/ref.chrom.sizes"
bedtools slop -i ${blacklist} -g ${chrsize} -b 500 | bedtools sort -i - | bedtools merge -i - > ATAC_bams/${tissue}/data/temp.bed
bedtools intersect -v -a ATAC_bams/${tissue}/idr_output/${mark}_${tissue}_overlap_peak.bed -b ATAC_bams/${tissue}/data/temp.bed | cut -f1-10 > ATAC_bams/${tissue}/data/peaks_no_blacklist.bed

#bedtools intersect -v -a ATAC_bams/${tissue}/idr_output/${mark}_${tissue}_overlap_peak.bed -b ${blacklist} | cut -f1-10 > ATAC_bams/${tissue}/data/peaks_no_blacklist.bed

#### Define train, validation and test chromosome splits
mkdir ATAC_bams/${tissue}/data/splits
module load conda 
conda activate chrombpnet
module load cuda/11.7.1 

module load seqtk
seqtk subseq refGenome/BosTaurusARSUCD12V105.chr.fa chrs.list | sed 's/>/>chr/g' > ATAC_bams/${tissue}/data/ref.fa
samtools faidx ATAC_bams/${tissue}/data/ref.fa
chrombpnet prep splits -c ${chrsize} -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ATAC_bams/${tissue}/data/splits/fold_0

#### Generate non-peaks (background regions)
rm -rf ATAC_bams/${tissue}/data/output_auxiliary
chrombpnet prep nonpeaks -g ATAC_bams/${tissue}/data/ref.fa -p ATAC_bams/${tissue}/data/peaks_no_blacklist.bed -c ${chrsize} -fl ATAC_bams/${tissue}/data/splits/fold_0.json -br ${blacklist} -o ATAC_bams/${tissue}/data/output
