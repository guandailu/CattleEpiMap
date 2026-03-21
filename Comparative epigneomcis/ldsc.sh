cd /group/zhougrp2/dguan/BovineFAANG/17_comparative_analysis/human_GWAS
mkdir -p Temp
##### make ann by RE type
for j in `cat tissue_specific_ann.list`
do
temp_ann_file=$(mktemp Temp/${j}.XXXXXX)
cat /group/zhougrp4/dguan/BovineFAANG/19_human_ref_CREs/05_tissue_specific/${j} | cut -f1-3 | awk '{print "chr"$0}' > ${temp_ann_file}
for i in {1..22}
do
python ldsc/make_annot.py --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim --bed-file ${temp_ann_file} --annot-file 03_annot/${j%.bed}.${i}.annot.gz
done
rm -rf ${temp_ann_file}
done
##### make ann by tissue
for j in `cat selected_tissues.txt`
do
temp_ann_file=$(mktemp Temp/${j}.XXXXXX)
cat /group/zhougrp4/dguan/BovineFAANG/19_human_ref_CREs/05_tissue_specific/${i}.E5.s*.tiss_spec.bed | cut -f1-3 | awk '{print "chr"$0}' > ${temp_ann_file}
for i in {1..22}
do
python ldsc/make_annot.py --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim --bed-file ${temp_ann_file} --annot-file 03_annot/${j}.${i}.annot.gz
done
done



##### calc l2 by RE type
for j in `cat tissue_specific_ann.list`
do
for i in {1..22}
do
python ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} --ld-wind-cm 1 --annot 03_annot/${j%.bed}.${i}.annot.gz --thin-annot --out 03_annot/${j%.bed}.${i} --print-snps hm3_no_MHC.list.txt
done
done
##### calc l2 by tiss
for j in `cat selected_tissues.txt`
do
for i in {1..22}
do
python ldsc/ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} --ld-wind-cm 1 --annot 03_annot/${j}.${i}.annot.gz --thin-annot --out 03_annot/${j}.${i} --print-snps hm3_no_MHC.list.txt
done
done



# format GWAS
### manually if no specific rules in GWAS summary statistic files

for i in AD Adhesions Afib AntThighFat Asthma BFP BMI BladderCa CD CRC CRF CerebVol Dementia DuodUlcer Estradiol GIBleed Gallstones GripStrength HCM HDL HF IBD Ileus KidneyCa KidneyVol LHipVol LiverAbscess LiverFe LiverVol LungAdCa LungCa LungVol MS MaleInfertility NonHDL PD PUD Pneumonia PostThighFat RA RHipVol SLE SpleenFe SpleenVol TC Testosterone UC
do
python ldsc/ldsc.py --h2 02_GWAS_reformatted/${i}.sumstats.gz --ref-ld-chr $(cat tissue_specific_ann.list | sed 's/.bed//g' | awk '{print "03_annot/"$0}' |  tr '\n' ' ' | sed 's/ /.,/g')baselineLD_v2.3/baselineLD. --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --overlap-annot --print-coefficients --print-delete-vals --out 05_h2_partitioned/${i}
done

