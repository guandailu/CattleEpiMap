for i in ../../13_trait_integration/01_GWAS_data/STgenetics.*.fdr.txt; 
do 
  t=$(echo $i | sed 's/.*STgenetics.//g;s/.fdr.txt//g')
  n=3530
  awk 'NR==FNR {a[$1]; next} FNR==1 || ($3 in a)' ref_snps_only.txt ../../13_trait_integration/01_GWAS_data/STgenetics.${t}.txt | awk '{if (NR == 1) print; else print "chr"$0}' > Temp/${t}.GWAS_filtered.txt
  python ldsc/munge_sumstats.py --sumstats Temp/${t}.GWAS_filtered.txt --N ${n} --out 01_gwas_data/${t} --a1 ref --a2 alt --snp SNP --p p 
done

module load bedtools2
for i in 02_annt_bed/*.bed
do
  a=$(echo $i | sed 's/.*\///g;s/.bed//g')
  echo "ANNOT" > 03_annt_make/${a}.annot
  bedtools intersect -a <(awk 'BEGIN{OFS="\t"} {print $1, $4-1, $4}' cattle_ref_panel.bim) -b 02_annt_bed/${a}.bed -c | awk '{if($4>0) print 1; else print 0}' >> 03_annt_make/${a}.annot
  gzip -f 03_annt_make/${a}.annot
done


for i in 03_annt_make/*.annot.gz
do
a=$(echo $i | sed 's/.*\///g;s/.annot.gz//g')
python ldsc/ldsc.py --l2 --bfile cattle_ref_panel --ld-wind-kb 1000 --annot 03_annt_make/${a}.annot.gz --thin-annot --out 03_annt_make/${a}
done
# for i in 03_annt_make/*.annot.gz; do sbatch -p bmh -c 4 --mem=48G -t 10-0 run_annot.sh $i; done


mkdir -p 04_ldsc_res
for i in 01_gwas_data/*.sumstats.gz
do
t=$(echo $i | sed 's/.*\///g;s/.sumstats.gz//g')
python ldsc/ldsc.py --h2 01_gwas_data/${t}.sumstats.gz --ref-ld 03_annt_make/genes,03_annt_make/exon,03_annt_make/intron,03_annt_make/five_prime_utr,03_annt_make/three_prime_utr,03_annt_make/conserved,03_annt_make/Cattle_eQTLs,03_annt_make/CattleQTLdb,03_annt_make/Kern_EnhA,03_annt_make/Kern_TssA,03_annt_make/Unique_EnhA_this,03_annt_make/Unique_TssA_this --w-ld 03_annt_make/genes --overlap-annot --frqfile cattle_ref_panel --print-coefficients --out 04_ldsc_res/${t}
done

cat 04_ldsc_res/Fat.results | awk '{if (NR==1) print $0"\ttrait"}' > 04_ldsc_res/Final_results_vsFAANG.sum.txt
for i in 04_ldsc_res/*.results; do t=$(echo $i | sed 's/.*\///g;s/.results//g'); cat $i | awk -v t=$t '{if (NR>1) print $0"\t"t}' >> 04_ldsc_res/Final_results_vsFAANG.sum.txt;done

mkdir -p 04_ldsc_res
for i in 01_gwas_data/*.sumstats.gz
do
t=$(echo $i | sed 's/.*\///g;s/.sumstats.gz//g')
python ldsc/ldsc.py --h2 01_gwas_data/${t}.sumstats.gz --ref-ld 03_annt_make/genes,03_annt_make/RM_top1pert --w-ld 03_annt_make/genes --overlap-annot --frqfile cattle_ref_panel --print-coefficients --out 04_ldsc_res/${t}
done



##### human traits
i="trait"
python ldsc/ldsc.py --h2 02_GWAS_reformatted/${i}.sumstats.gz --ref-ld-chr $(cat tissue_specific_ann.list | sed 's/.bed//g' | awk '{print "03_annot/"$0}' |  tr '\n' ' ' | sed 's/ /.,/g')baselineLD_v2.3/baselineLD. --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --overlap-annot --print-coefficients --print-delete-vals --out 05_h2_partitioned/${i}

