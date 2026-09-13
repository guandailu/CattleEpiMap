ldetect run --genetic-map ldetect/out/chr${i}_map.gz --reference-panel ldetect/out/chr${i}.vcf.gz --individuals ldetect/out/Holstein_indiv.list --chromosome ${i} --output-dir results/chr${i}
#!/bin/bash

i=$1
mkdir -p 01_GWAS_data1
#for i in ../13_trait_integration/01_GWAS_data/USdairy.*.fdr.txt;
#do
#i="../13_trait_integration/01_GWAS_data/USdairy.Teat_length.txt"
t=$(echo $i | sed 's/.*USdairy.//g;s/.fdr.txt//g')
mkdir -p 02_matched_ldm/${t}
n=$(cat /group/zhougrp2/dguan/BovineFAANG/revision/01_LDSC/trait_summary.csv  | awk -v t=$t -F "," '{if ($2 == t) print $3}')

if [ -z "$n" ]; then
    echo "ERROR: Could not find sample size for trait: $t"
else
    echo "Processing $t with N=$n..."
    
    awk -v n="$n" '
      # 1. Process the MAF file first
      NR==FNR {
        if(NR>1) {
          maf[$1] = $4  # $1 is SNP ID (e.g. 1:6825), $4 is MAF
        }
        next
      }
      
      # 2. Print the GCTB Header
      FNR==1 { 
        print "SNP", "A1", "A2", "MAF", "b", "se", "p", "N" 
      }
      
      # 3. Process the GWAS file and match by SNP ID
      FNR>1 {
        # GWAS file columns: $3=SNP, $4=ref, $5=alt, $6=p, $7=beta, $8=se
        if ($3 in maf) {
          # GCTB format: SNP A1 A2 MAF b se p N
          # We use $5 (alt) as A1 and $4 (ref) as A2 as is standard for GCTB
          print $3, $5, $4, maf[$3], $7, $8, $6, n
        }
      }
    ' OFS="\t" 12M.af.info.tsv ../13_trait_integration/01_GWAS_data/USdairy.${t}.txt > 01_GWAS_data1/USDairy.${t}.gctb.txt
fi
#done

./gctb_2.5.5_Linux/gctb --bfile refpanel.filtered --make-block-ldm --block-info results/Holstein_ld_blocks_fixed.bed --out ldm --thread 24

./gctb_2.5.5_Linux/gctb --ldm ldm --make-ldm-eigen --out ldm --thread 12

./gctb_2.5.5_Linux/gctb --thread 1 --ldm ldm --gwas-summary 01_GWAS_data/USDairy.${t}.gctb.txt --make-ldm-eigen --out 02_matched_ldm/${t}/

t="Body_depth"
./gctb_2.5.5_Linux/gctb --gwfm RC --ldm-eigen 02_matched_ldm/${t}/ --gwas-summary 01_GWAS_data/USDairy.${t}.gctb.txt --annot 03_annt/RM_calib.txt --thread 1 --out 04_fine_mapped/${t}

