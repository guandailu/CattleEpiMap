library(data.table)
library(tidyverse)
library(rstatix)
ase = fread("/group/zhougrp2/Wenjing/0_AS_run/4_AS/Summary_AS_results/AS_0.05_M08_22.bed", header=T) %>% mutate(SNP=paste(chr, position, ref, alt, sep="_"), logfc = log2(altCount / refCount)) %>% separate(sample_id, into=c("Assay", "Tissue", "ID"), sep="_", remove=F) %>% filter(qvalue < 0.05) %>% select(SNP, Assay, Tissue, logfc) 
assays = c("ATAC", "CTCF", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")
svm=data.frame()
for (a in assays){
  tissues=fread(paste0("/group/zhougrp4/dguan/BovineFAANG/23_deep_learning/lsgkm/", a,"_tissues.list"), header=F) %>% pull(V1) %>% as.vector
  for (t in tissues){
    df = fread(paste0("ASEout/", a,"/", t,".ASE_tested_SNPs.deltaSVM.tsv"))
    names(df)=c("SNP", "deltaSVM")
    df %>% mutate(Assay = a, Tissue = t) %>% select(SNP, Assay, Tissue, deltaSVM) -> df
    svm = rbind(svm, df)
    rm(df)
  }
}

saveRDS(svm, file="ASEout/ASE_tested_SNPs.deltaSVM.RDS")
merge(ase, svm, by=c("SNP", "Assay", "Tissue")) %>% distinct -> mydf
svm %>% group_by(Tissue) %>% summarize(num_assays=length(unique(Assay))) %>% filter(num_assays == 6) %>% pull(Tissue) -> tissues
for (a in assays){
  for (t in tissues){
    subset(mydf, Assay == a & Tissue == t) -> plot_df
    if (nrow(plot_df)> 50){
      p=ggplot(plot_df)+
            geom_point(aes(x=logfc, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=logfc, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("log2FC") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
      ggsave(paste0("ASEout/", a, ".", t, ".log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)
      plot_df %>% cor_test(logfc, deltaSVM) %>% mutate(Assay=a, Tissue = t) -> cor_df
      fwrite(cor_df, paste0("ASEout/", a, ".", t, ".log2FC_vs_deltaSVM.txt"), sep="\t", row.names=F, quote=F)
    }
  }
}



library(data.table)
library(tidyverse)
library(rstatix)
ase = fread("AS_0.05_M08_22.bed", header=T) %>% mutate(SNP=paste(chr, position, ref, alt, sep="_"), logfc = log2(altCount / refCount)) %>% separate(sample_id, into=c("Assay", "Tissue", "ID"), sep="_", remove=F) %>% filter(qvalue < 0.05) %>% select(SNP, Assay, Tissue, logfc) 
assays = c("ATAC", "CTCF", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")
svm=readRDS("ASEout/ASE_tested_SNPs.deltaSVM.RDS")
merge(ase, svm, by=c("SNP", "Assay", "Tissue")) %>% distinct -> mydf
a="ATAC"
tissues=fread("../chrombpnet/final_tissues_with_chrombpnet_models.txt", header=F) %>% pull(V1) %>% as.vector
for (t in tissues){
subset(mydf, Assay == a & Tissue == t) -> sub_df
chrombpnet = fread("../chrombpnet/ATAC_bams/Kidney/ASE_scores/ASE_tested_SNPs.list_snp_scores.tsv", header=T) %>% mutate(SNP=paste(str_replace(CHR, "chr",""), POS0, REF, ALT, sep="_")) 
merge(sub_df, chrombpnet, by = "SNP") -> plot_df
if (nrow(plot_df) > 10){
print(t)
#plot_df %>% mutate(abs_log_counts_diff = abs(log_counts_diff), abs_deltaSVM = abs(deltaSVM)) %>% cor_test(abs_log_counts_diff, abs_deltaSVM) %>% print

plot_df %>% mutate(log_counts_diff = ifelse(log_counts_diff > 0, log_counts_diff * -1, abs(log_counts_diff))) %>% cor_test(log_counts_diff, deltaSVM) %>% print
}
}

a="ATAC"
t="Rumen"
subset(mydf, Assay == a & Tissue == t) -> sub_df
chrombpnet = fread(paste0("../chrombpnet/ATAC_bams/", t,"/ASE_scores/ASE_tested_SNPs.list_snp_scores.tsv"), header=T) %>% mutate(SNP=paste(str_replace(CHR, "chr",""), POS0, REF, ALT, sep="_")) 
merge(sub_df, chrombpnet, by = "SNP") -> plot_df
if (nrow(plot_df) > 10){
	print(t)
	#plot_df %>% mutate(abs_log_counts_diff = abs(log_counts_diff), abs_deltaSVM = abs(deltaSVM)) %>% cor_test(abs_log_counts_diff, abs_deltaSVM) %>% print
  #plot_df %>% mutate(log_counts_diff = ifelse(log_counts_diff > 0, log_counts_diff * -1, abs(log_counts_diff))) %>% cor_test(log_counts_diff, deltaSVM) %>% print
  plot_df %>%  cor_test(log_counts_diff, deltaSVM) %>% print
}
p=ggplot(plot_df)+
            geom_point(aes(x=log_counts_diff, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=log_counts_diff, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("ChromBPNet effect size\n(log2FC)") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave(paste0("ASEout/", a, ".", t, ".ChromBPNet_log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)



