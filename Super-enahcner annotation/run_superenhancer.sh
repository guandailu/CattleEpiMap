#!/bin/bash
# this is by individual epigenome, will not consider in future analysis
#run_superenhancer.sh
set -e

### parse parameters
eid=$1
path="./"
### load mandatory modules
module load samtools
module load R
module load bedtools2

### create temporary directory and file
mkdir -p ${path}/04_super_enhacner/ROSE/Temp
mkdir -p ${path}/04_super_enhacner/output/${eid}
enhancer_file=$(mktemp ${path}/04_super_enhacner/ROSE/Temp/enhancer.XXXXXXXXXXXXXX.gff)

### reheader bam file
samtools view -h ${path}/01_data_preprocess/05_mappability_filter/H3K27ac_${eid}.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam
samtools index ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam

### prepare enhancer file
for e in E{4..6}
do
    cat ${path}/03_chromatin_ann/FourMarks/ChromatinStates/${eid}_States.bed | awk -v e=$e '{if ($4 == e) print}' | bedtools merge -i - | awk -v OFS="\t" '{print "chr"$1,$1"_"$2"_"$3,".",$2,$3,".", ".", ".", $1"_"$2"_"$3}' >> ${enhancer_file}
done

### run ROSS
PATHTO="${path}/04_super_enhacner/ROSE"
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin
python3 ${PATHTO}/bin/ROSE_main.py -g ARSUCD12 -i ${enhancer_file} -r ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam -o ${path}/04_super_enhacner/output/${eid}

### remove temporary files
wait
rm -rf ${path}/04_super_enhacner/ROSE/Temp/H3K27ac_${eid}.renameChr.bam ${enhancer_file}
echo -e "Done!!!\n"

# for i in `cat ../../02_imputation/imputed_samples.tab | cut -f1 | sort | uniq`; do sbatch -p high -c 4 --mem=12G -t 10-0 -o Logs/superenh.${i}.%j.out run_superenhancer.sh $i; done
