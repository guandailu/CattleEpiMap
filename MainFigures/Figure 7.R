library(data.table)
library(tidyverse)
library(rstatix)

files=list.files("../../cattleGTEx/Zng3-3ElohU/", pattern=".log2aFC.txt.gz$")
eqtl_df = data.frame()
for (f in 1:length(files)){
  df = fread(paste0("../../cattleGTEx/Zng3-3ElohU/", files[f]), header=T) %>% mutate(Tissue=str_replace(files[f], ".log2aFC.txt.gz","")) %>% select(Tissue, pheno_id, variant_id, log2_aFC, start_distance, af) 
  eqtl_df = rbind(eqtl_df, df)
}

assays = c("ATAC", "CTCF", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")
svm=data.frame()
for (a in assays){
  tissues=fread(paste0("/group/zhougrp4/dguan/BovineFAANG/23_deep_learning/lsgkm/", a,"_tissues.list"), header=F) %>% pull(V1) %>% as.vector
  for (t in tissues){
    df = fread(paste0("eQTLout/", a,"/", t,".eQTL_SNP.deltaSVM.tsv"))
    names(df)=c("SNP", "deltaSVM")
    df %>% mutate(Assay = a, Tissue = t) %>% select(SNP, Assay, Tissue, deltaSVM) -> df
    svm = rbind(svm, df)
    rm(df)
  }
}
saveRDS(svm, file="eQTLout/eQTL_SNP.deltaSVM.RDS")
svm %>% separate(SNP, into=c("chr", "pos", "ref", "alt"), remove=F) %>% mutate(SNP=paste0(chr, "_", pos)) %>% select(SNP, Assay, Tissue, deltaSVM) -> svm
names(eqtl_df)[3]="SNP"
merge(eqtl_df, svm, by = c("SNP", "Tissue")) -> mydf
tissues=intersect(unique(eqtl_df$Tissue), unique(svm$Tissue))

mydf %>% filter(Assay == "ATAC") %>% filter(log2_aFC > 3 | log2_aFC < -3) %>% cor_test(log2_aFC, deltaSVM) %>% print

library(rstatix)

a="ATAC"
for (t in tissues){
	subset(mydf, Assay == a & Tissue == t) -> plot_df
	if (nrow(plot_df) > 10){
		print(t)
  	print(nrow(plot_df))
		#plot_df %>% mutate(abs_log_counts_diff = abs(log_counts_diff), abs_deltaSVM = abs(deltaSVM)) %>% cor_test(abs_log_counts_diff, abs_deltaSVM) %>% print
		plot_df %>% filter(log2_aFC > 1 | log2_aFC < -1) %>% cor_test(log2_aFC, deltaSVM) %>% print
	}
}

a="ATAC"
t="Epididymis"
subset(mydf, Assay == a & Tissue == t) %>% filter(log2_aFC > 1 | log2_aFC < -1) %>% filter(log2_aFC > -6.643856) -> plot_df
plot_df %>% cor_test(log2_aFC, deltaSVM) %>% print
p=ggplot(plot_df)+
            geom_point(aes(x=log2_aFC, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=log2_aFC, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("log2_aFC") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave(paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)
plot_df %>% cor_test(log2_aFC, deltaSVM) -> cor_df
fwrite(cor_df, paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.txt"), sep="\t", row.names=F, quote=F)



a="ATAC"
t="Adipose"
subset(mydf, Assay == a & Tissue == t) %>% filter(log2_aFC > 1 | log2_aFC < -1)  -> plot_df
plot_df %>% cor_test(log2_aFC, deltaSVM) %>% print
p=ggplot(plot_df)+
            geom_point(aes(x=log2_aFC, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=log2_aFC, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("log2_aFC") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave(paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)
plot_df %>% cor_test(log2_aFC, deltaSVM) -> cor_df
fwrite(cor_df, paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.txt"), sep="\t", row.names=F, quote=F)



a="ATAC"
t="Liver"
subset(mydf, Assay == a & Tissue == t) %>% filter(log2_aFC > 1 | log2_aFC < -1)  -> plot_df
plot_df %>% cor_test(log2_aFC, deltaSVM) %>% print
p=ggplot(plot_df)+
            geom_point(aes(x=log2_aFC, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=log2_aFC, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("log2_aFC") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave(paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)
plot_df %>% cor_test(log2_aFC, deltaSVM) -> cor_df
fwrite(cor_df, paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.txt"), sep="\t", row.names=F, quote=F)



a="ATAC"
t="Muscle"
subset(mydf, Assay == a & Tissue == t) %>% filter(log2_aFC > 1 | log2_aFC < -1)  -> plot_df
plot_df %>% cor_test(log2_aFC, deltaSVM) %>% print
p=ggplot(plot_df)+
            geom_point(aes(x=log2_aFC, y = deltaSVM), size=2, alpha=0.6, color="black")+
            geom_smooth(aes(x=log2_aFC, y = deltaSVM), method="lm", color="red", se=F)+
            xlab("log2_aFC") +
            ylab("deltaSVM")+
            theme_classic(base_size = 20, base_line_size = 0.5)+
            theme(
            axis.text.x = element_text(color="black", size=15),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave(paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.pdf"), p, width = 3.5, height=3)
plot_df %>% cor_test(log2_aFC, deltaSVM) -> cor_df
fwrite(cor_df, paste0("eQTLout/", a, ".", t, ".eQTL_log2FC_vs_deltaSVM.txt"), sep="\t", row.names=F, quote=F)
