#!/bin/bash

set -e
#run_validation.sh

batch=$1

mkdir -p cs${batch} cs${batch}/Temp/ cs${batch}/Logs/
 
### convert data format
jobids=()
num=$(cat imputed_samples.sample_with_seed${batch}.tab | wc -l)
for i in `seq $num`
do 
  cat imputed_samples.sample_with_seed${batch}.tab | awk -v i=$i '{if (NR == i) print}' > cs${batch}/Temp/sample_line${i}.tab
  job=$(sbatch -p bmm -c 1 --mem=4G -t 10-0 --job-name=cs${batch}.Convert -o cs${batch}/Logs/convert.line${i}.%j.out Convert.sh cs${batch}/Temp/sample_line${i}.tab ${batch} | awk '{print $4}')
  jobids+=(${job})
done
all_jobs=$(echo ${jobids[*]} | sed 's/ /:/g')


### ComputeGlobalDist
jobids=()
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF
do
  job=$(sbatch -p bmh -c 8 --mem=32G -t 10-0 --job-name=cs${batch}.ComputeGlobalDist --dependency=afterany:${all_jobs} -o cs${batch}/Logs/ComputeGlobalDist.${i}.%j.out ComputeGlobalDist.sh $i imputed_samples.sample_with_seed${batch}.tab ${batch} | awk '{print $4}')
  jobids+=(${job})
done
all_jobs=$(echo ${jobids[*]} | sed 's/ /:/g')

jobids=()
for i in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF
  do 
  jobn=$(sbatch -p bmh -c 8 --mem=32G -t 10-0 --job-name=cs${batch}.GenerateTrainData --dependency=afterany:${all_jobs} -o cs${batch}/Logs/GenerateTrainData.${i}.%j.out GenerateTrainData.sh $i imputed_samples.sample_with_seed${batch}.tab ${batch} | awk '{print $4}')
  jobids+=(${jobn})
done
all_jobs=$(echo ${jobids[*]} | sed 's/ /:/g')


for i in `cat imputed_samples.sample_with_seed${batch}.tab | awk '{print $1}' | sort | uniq`
  do 
  for j in H3K4me1 H3K4me3 H3K27ac H3K27me3 CTCF
    do sbatch -p med -c 4 --mem=8G -t 10-0 --job-name=cs${batch}.TrainApply --dependency=afterany:${all_jobs} -o cs${batch}/Logs/GeneratePredictors.${i}_${j}.%j.out TrainApply.sh $i $j imputed_samples.sample_with_seed${batch}.tab ${batch}
  done
done
