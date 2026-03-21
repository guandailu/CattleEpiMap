#!/bin/bash
set -e

sample_tab=$1
module load openjdk/16.0.2
mkdir -p CONVERTEDDIR
java -jar -Xmx4G /home/dguan/bin/ChromImpute/ChromImpute.jar Convert /group/zhougrp4/dguan/BovineFAANG/02_imputation/01_signal_tracks/ ${sample_tab} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome CONVERTEDDIR
#num=$(cat imputed_samples.tab | wc -l); for i in `seq $num`; do cat imputed_samples.tab | awk -v i=$i '{if (NR == i) print}' > Temp/sample_line${i}.tab; sbatch -p bmm -c 1 --mem=4G -t 10-0 -o Logs/convert.line${i}.%j.out ChromImpute_Convert.sh Temp/sample_line${i}.tab; done

#!/bin/bash
set -e
module load openjdk/16.0.2
mark=$1
sample_tab=$2
mkdir -p DISTANCEDIR
java -jar -Xmx32G /home/dguan/bin/ChromImpute/ChromImpute.jar ComputeGlobalDist -m ${mark} CONVERTEDDIR ${sample_tab} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome DISTANCEDIR
java -jar -Xmx32G /home/dguan/bin/ChromImpute/ChromImpute.jar GenerateTrainData CONVERTEDDIR DISTANCEDIR ${sample_tab} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome TRAINDATA ${mark}
## Compute global distance
#for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF; do sbatch -p bmh -c 8 --mem=32G -t 10-0 ChromImpute_ComputeGlobalDist.sh $i imputed_samples.tab; done

#!/bin/bash
set -e
module load openjdk/16.0.2
EID=$1
assay=$2
sample_tab=$3
java -jar -Xmx8G /home/dguan/bin/ChromImpute/ChromImpute.jar Train TRAINDATA ${sample_tab} PREDICTORDIR ${EID} ${assay}
java -jar -Xmx8G /home/dguan/bin/ChromImpute/ChromImpute.jar Apply CONVERTEDDIR DISTANCEDIR PREDICTORDIR ${sample_tab} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome IMPUTED ${EID} ${assay}
#Generate the trained predictors
#for i in `cat imputed_samples.tab | awk '{print $1}' | sort | uniq`; do for j in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF; do sbatch -p high -c 4 --mem=8G -t 10-0 -o Logs/GeneratePredictors.${i}_${j}.%j.out ChromImpute_TrainApply.sh $i $j imputed_samples.tab; done; done

#!/bin/bash
set -e
module load openjdk/16.0.2
sample_tab=$1
mkdir -p CHROMHMMDIR
java -jar -Xmx24G /home/dguan/bin/ChromImpute/ChromImpute.jar ExportToChromHMM -g 2 IMPUTED ${sample_tab} /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome CHROMHMMDIR
# Export to chromHMM
#for i in `cat imputed_samples.tab | awk '{print $1}' | sort | uniq`; do for j in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF; do echo -e "${i}\t${j}\timpute_${i}_${j}.wig.gz" >> export2chromhmm.tab; done; done
#num=$(cat export2chromhmm.tab | wc -l); for i in `seq $num`; do cat export2chromhmm.tab | awk -v i=$i '{if (NR == i) print}' > Temp/export2chromhmm.${i}.tab; sbatch -p bmm -c 2 --mem=6G -t 10-0 -o Logs/export2chromhmm.${i}.%j.out Export2ChromHMM.sh Temp/export2chromhmm.${i}.tab; done
