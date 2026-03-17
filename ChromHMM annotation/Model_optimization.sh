mkdir -p Rep1/BINARYDIR Rep2/BINARYDIR
cp -r BINARYDIR/*M08* Rep1/BINARYDIR/.
cp -r BINARYDIR/*M22* Rep2/BINARYDIR/.

### we run model from 2 to 20
#!/bin/bash
# Script: ChromHMM_ModelOptm.sh

module load openjdk/16.0.2
mkdir -p Logs

num_model=$1
java -jar -Xmx64G /home/dguan/bin/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome Rep1/BINARYDIR Rep1/LearnModel_${num_model} ${num_model} bosTau9
java -jar -Xmx64G /home/dguan/bin/ChromHMM/ChromHMM.jar LearnModel -p 12 -l /group/zhougrp3/dguan/BovineFAANG/refGenome/BosTaurusARSUCD12V105.chr.fa.genome Rep2/BINARYDIR Rep1/LearnModel_${num_model} ${num_model} bosTau9
# for i in {2..20}; do sbatch -p bmm -c 12 --mem=48G -t 10-0 -o LearnModel_${i}.%j.out ChromHMM_LearnModel.sh $i; done
