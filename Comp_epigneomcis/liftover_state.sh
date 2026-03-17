
HumanStateFolder="17_comparative_analysis/01_human_state_bytissue/"
CattleStateFolder="17_comparative_analysis/02_cattle_state_bytissue"

state=$1
cattle_state_file="${CattleStateFolder}/${state}_consensus.bed"
tmp_file1=$(mktemp Temp/cattle_CRE.XXXXXXXXX.bed)
cat $cattle_state_file | awk '{print "chr"$0,$1"_"$2"_"$3}' > ${tmp_file1}
lift_output=$(mktemp Temp/lift_CRE.XXXXXXXXX.bed)
unlift_output=$(mktemp Temp/unlift_CRE.XXXXXXXXX.bed)
### lift cattle CRE to human genome
/home/dguan/bin/ucsc_utilities/liftOver -minMatch=0.1 -bedPlus=5 ${tmp_file1} bosTau9ToHg19.over.chain ${lift_output} ${unlift_output}

#### cattle species specific CRE
cat ${unlift_output} | grep -v "^#" | bedtools sort -i - | sed 's/chr//g' | awk -v OFS="\t" 'BEGIN{print "CattleCREid\tCREchr\tCREstart\tCREend"}{print $4,$1,$2,$3}' > 03_cattle2human_CREs/${state}.cattle2human.ssCRE.bed


### sequence and functional conserved CREs
human_state_file="${HumanStateFolder}/${state}_consensus.bed"
bedtools intersect -a <(cat ${lift_output} | sed 's/chr//g' | bedtools sort -i -) -b <(cat ${human_state_file} | bedtools sort -i - ) -f 0.5 -wa -wb | awk -v OFS="\t" 'BEGIN{print "CattleCREid\tHumanCREchr\tHumanCREstart\tHumanCREend\tMatchHumanChr\tMatchHumanStart\tMatchHumanEnd"}{print $4,$5,$6,$7,$1,$2,$3}' > 03_cattle2human_CREs/${state}.cattle2human.sfCRE.bed
seqcons_output=$(mktemp Temp/sc_CRE.XXXXXXXXX.bed)
grep -v -f <(cat 03_cattle2human_CREs/${state}.cattle2human.sfCRE.bed | awk '{if (NR >= 2) print $1}' | sort) <(cat ${lift_output} | sort -k4,4) > ${seqcons_output}
    
### sequence conserved but functional divergent CREs
elements=("E1" "E2" "E3" "E4" "E5" "E6" "E7" "E8" "E9" "E10" "E11" "E12" "E13" "E14")
rest_files=( $(echo ${elements[@]} | tr ' ' '\n' | grep -vEw "$state" | awk '{print "01_human_state_bytissue/"$1"_consensus.bed"}' | tr '\n' ' ') )
bedtools intersect -a <(cat ${seqcons_output} | sed 's/chr//g' | bedtools sort -i -) -b <(cat ${rest_files} | bedtools sort -i -) -f 0.5 -wa -wb | awk -v OFS="\t" 'BEGIN{print "CattleCREid\tHumanCREchr\tHumanCREstart\tHumanCREend\tMatchHumanChr\tMatchHumanStart\tMatchHumanEnd"}{print $4,$5,$6,$7,$1,$2,$3}' | uniq > 03_cattle2human_CREs/${state}.cattle2human.sdCRE.bed

### only sequence conserved CREs
grep -v -f <(cat 03_cattle2human_CREs/${state}.cattle2human.sdCRE.bed | awk '{if (NR >= 2) print $1}' | sort) <(cat ${seqcons_output} | sed 's/chr//g' | bedtools sort -i -) | awk -v OFS="\t" 'BEGIN{print "CattleCREid\tMatchHumanChr\tMatchHumanStart\tMatchHumanEnd"}{print $4,$1,$2,$3}'  > 03_cattle2human_CREs/${state}.cattle2human.soCRE.bed
rm -rf ${tmp_file1} ${lift_output} ${unlift_output} ${seqcons_output}


