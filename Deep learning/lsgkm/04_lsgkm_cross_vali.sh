#!/bin/bash

mark=$1
tissue=$2
cv=$3

#cd ${tissue}
mkdir -p output/${mark}/${tissue}/cross_vali

gkmtrain -x 10 -i ${cv} -T 16 output/${mark}/${tissue}/${mark}_${tissue}_overlap_peak.pos1x.fa output/${mark}/${tissue}/${mark}_${tissue}_overlap_peak.neg1x.fa output/${mark}/${tissue}/cross_vali/${mark}_${tissue}.${cv}
