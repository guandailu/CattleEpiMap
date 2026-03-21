#!/bin/bash
#run_superenhancer.byTissue.sh
set -e

### parse parameters
tissue=$1
path="./"
### load mandatory modules
module load samtools
module load R
module load bedtools2

### create temporary directory and file
mkdir -p ${path}/04_super_enhacner/ROSE/Temp
mkdir -p ${path}/04_super_enhacner/output/${tissue}
enhancer_file=$(mktemp ${path}/04_super_enhacner/ROSE/Temp/enhancer.XXXXXXXXXXXXXX.gff)

### get eid for a tissue
eids=( $(cat ${path}/02_imputation/imputed_samples.tab | grep -E "^${tissue}" | cut -f1 | sort | uniq | tr '\n' ' ') )
echo "Epigenomes for ${tissue} include: ${eids[@]}"

### reheader bam file
for eid in ${eids[@]}
    do
    samtools view -h ${path}/01_data_preprocess/05_mappability_filter/H3K27ac_${eid}.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam
    samtools index ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam

    samtools view -h ${path}/01_data_preprocess/05_mappability_filter/Input_${eid}.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > ${path}/04_super_enhacner/ROSE/Temp/Input_${eid}.renameChr.bam
    samtools index ${path}/04_super_enhacner/ROSE/Temp/Input_${eid}.renameChr.bam
done

### bam by tissue
mkdir -p BamByTissue
bams=$(cat ${path}/02_imputation/imputed_samples.tab | grep -E "^${tissue}" | cut -f1 | sort | uniq | awk -v p=${path} '{print p"/04_super_enhacner/ROSE/Temp/H3K27ac_"$1".renameChr.bam"}' | tr '\n' ' ')
samtools merge -f -@ 8 -l 9 -o ${path}/04_super_enhacner/ROSE/BamByTissue/H3K27ac_${tissue}.bam ${bams}
samtools index ${path}/04_super_enhacner/ROSE/BamByTissue/H3K27ac_${tissue}.bam
inputs=$(echo $bams | sed 's/H3K27ac_/Input_/g')
samtools merge -f -@ 8 -l 9 -o ${path}/04_super_enhacner/ROSE/BamByTissue/Input_${tissue}.bam ${inputs}
samtools index ${path}/04_super_enhacner/ROSE/BamByTissue/Input_${tissue}.bam

### prepare enhancer file
cat `cat ${path}/02_imputation/imputed_samples.tab | grep -E "^${tissue}" | cut -f1 | sort | uniq | awk -v p=${path} '{print p"/03_chromatin_ann/FourMarks/ChromatinStates/"$1"_States.bed"}' | tr '\n' ' '` | awk '{if ($4 == "E4" || $4 == "E5" || $4 == "E6") print}' | bedtools sort -i - | bedtools merge -i - | awk -v OFS="\t" '{print "chr"$1,$1"_"$2"_"$3,".",$2,$3,".", ".", ".", $1"_"$2"_"$3}' > ${enhancer_file}
# here 250 = 534.199 / 2, 534.199 is the avergae length of enhancers
### run ROSS
PATHTO="${path}/04_super_enhacner/ROSE"
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
mkdir -p ${path}/04_super_enhacner/output/${tissue}
python3 ${PATHTO}/bin/ROSE_main.py -g ARSUCD12 -i ${enhancer_file} -r ${path}/04_super_enhacner/ROSE/BamByTissue/H3K27ac_${tissue}.bam -c ${path}/04_super_enhacner/ROSE/BamByTissue/Input_${tissue}.bam -o ${path}/04_super_enhacner/output/${tissue}

### remove temporary files
wait
echo -e "Done!!!\n"
