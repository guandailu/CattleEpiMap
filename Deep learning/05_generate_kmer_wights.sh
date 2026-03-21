#!/bin/bash

#python scripts/nrkmers.py 11 data/11_mer.fa

#./src/gkmpredict data/11_mer.fa data/ATAC/ATAC_lsgkm_aggregated_models.txt data/11_mer.weights.txt


#for i in `cat ../chrombpnet/ATAC_tissues.list`
#do
# echo -e '#!/bin/bash\n' "./src/gkmpredict data/11_mer.fa data/ATAC/${i}/ATAC_${i}.model.txt data/ATAC/${i}/ATAC_${i}.weights.txt" | sbatch -p med -c 2 --mem=6G -t 10-0 -o Logs/gkm.weight.${i}.%j.out
#done


for a in ATAC CTCF H3K27ac H3K27me3 H3K4me1 H3K4me3
do
    for i in `cat ${a}_tissues.list`
    do
        echo -e '#!/bin/bash\n' "./src/gkmpredict data/11_mer.fa output/${a}/${i}/${a}_${i}.model.txt output/${a}/${i}/${a}_${i}.weights.txt" | sbatch -p bmm -c 2 --mem=6G -t 3-0 -o Logs/gkm.weight.${a}.${i}.%j.out
    done
done
