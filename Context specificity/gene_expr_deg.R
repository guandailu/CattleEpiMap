library(data.table)
library(tidyverse)
library(DESeq2)
library("apeglm")

#### load data
df1 = fread("featureCounts_SE_output.tsv") %>% as.data.frame 
df1 %>% select(1,7:ncol(df1)) -> df1
df2 = fread("featureCounts_output.tsv") %>% as.data.frame 
df2 %>% select(1,7:ncol(df2)) -> df2
merge(df1, df2, by = "Geneid") -> df

#### prepare meta data
meta1=data.frame(Sample=as.data.frame(colnames(df1)[2:ncol(df1)]), seq_end="single")
names(meta1)=c("BioSample", "SeqEnd")
meta2=data.frame(Sample=as.data.frame(colnames(df2)[2:ncol(df2)]), seq_end="paired")
names(meta2)=c("BioSample", "SeqEnd")
meta_tab = rbind(meta1, meta2)


#### wroking for muscle
cts_df = df[,c("Geneid", "FetalMuscle_410", "FetalMuscle_438", "FetalMuscle_500", "FetalMuscle_503", "FetalSemimembranosusMuscle_2181", "FetalSemimembranosusMuscle_6819", "Muscle_Daisy", "Muscle_M08", "Muscle_M22", "SemimembranosusMuscle_6819", "SemimembranosusMuscle_2181", "LongissimusThoracisMuscle_6819", "LongissimusThoracisMuscle_2181")]
# prepare data for DESeq2
colnames(cts_df)[2:ncol(cts_df)] %>% as.data.frame -> samples
names(samples)="BioSample"
merge(samples, meta_tab, by= "BioSample") -> samples
samples$source = ifelse(grepl("_M",samples$BioSample),"US","AUS")
samples$condition=factor(ifelse(grepl("Fetal", samples$BioSample), "Fetal", "Adult"))
rownames(samples)=samples$BioSample
samples$BioSample=NULL
rownames(cts_df)=cts_df$Geneid
cts_df$Geneid=NULL
# run DESeq2
cts_df=cts_df[ ,order(names(cts_df))]
samples = samples[order(rownames(samples)),]
dds <- DESeqDataSetFromMatrix(countData = as.matrix(cts_df),
                              colData = as.matrix(samples),
                              design = ~ source + condition)
keep <- rowSums(counts(dds) >= 10) >= 3 #at least 3 samples with a count of 10 or more
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, name="condition_Fetal_vs_Adult")
res <- results(dds, contrast=c("condition","Fetal","Adult"))
resLFC <- lfcShrink(dds, coef="condition_Fetal_vs_Adult", type="apeglm")
resOrdered <- resLFC[order(resLFC$pvalue),]
fwrite(as.data.frame(resOrdered), "Muscle.fetal_vs_adult.DEGs.txt", sep="\t", row.names=T, quote=F)


          
myDEA=function(tissue, cts_df){
    # prepare data for DESeq2
    colnames(cts_df)[2:ncol(cts_df)] %>% as.data.frame -> samples
    names(samples)="BioSample"
    merge(samples, meta_tab, by= "BioSample") -> samples
    samples$source = ifelse(grepl("6819|2181|Daisy",samples$BioSample),"AUS","US")
    samples$condition=factor(ifelse(grepl("Fetal", samples$BioSample), "Fetal", "Adult"))
    subset(samples, condition == "Fetal") %>% pull(BioSample) -> fetalSamples
    cat("Fetal samples: ", fetalSamples, "\n")
    subset(samples, condition == "Adult") %>% pull(BioSample) -> adultSamples
    cat("Adult samples: ", adultSamples, "\n")
    rownames(samples)=samples$BioSample
    samples$BioSample=NULL
    rownames(cts_df)=cts_df$Geneid
    cts_df$Geneid=NULL
    # run DESeq2
    cts_df=cts_df[ ,order(names(cts_df))]
    samples = samples[order(rownames(samples)),]
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(cts_df),
                                  colData = as.matrix(samples),
                                  design = ~ source + condition)
    keep <- rowSums(counts(dds) >= 10) >= 3 #at least 3 samples with a count of 10 or more
    dds <- dds[keep,]
    dds <- DESeq(dds)
    res <- results(dds)
    res <- results(dds, name="condition_Fetal_vs_Adult")
    res <- results(dds, contrast=c("condition","Fetal","Adult"))
    resLFC <- lfcShrink(dds, coef="condition_Fetal_vs_Adult", type="apeglm")
    resOrdered <- resLFC[order(resLFC$pvalue),]
    fwrite(as.data.frame(resOrdered), paste0(tissue, ".fetal_vs_adult.DEGs.txt"), sep="\t", row.names=T, quote=F)
}


#### wroking for Liver
cts_df = df[,c("Geneid", "FetalLiver_410", "FetalLiver_438", "FetalLiver_500", "FetalLiver_503", "FetalLiver_2181", "FetalLiver_6819", "Liver_2181", "Liver_6819", "Liver_Daisy", "Liver_M08", "Liver_M22")]
myDEA("Liver", cts_df)

#### wroking for Heart
colnames(df)[grepl("Fetal", colnames(df))]
colnames(df)[grepl("Heart", colnames(df))]
cts_df = df[,c("Geneid", "FetalHeart_410", "FetalHeart_500", "FetalHeart_503", "FetalHeart_2181", "FetalHeart_6819", "Heart_2181", "Heart_6819", "Heart_Daisy")]
myDEA("Heart", cts_df)

#### wroking for Lung
colnames(df)[grepl("Fetal", colnames(df))]
colnames(df)[grepl("Lung", colnames(df))]
cts_df = df[,c("Geneid", "FetalLung_2181", "FetalLung_6819", "Lung_2181", "Lung_6819", "Lung_Daisy", "Lung_M08", "Lung_M22")]
myDEA("Lung", cts_df)

#### wroking for Kidney
colnames(df)[grepl("Fetal", colnames(df))]
colnames(df)[grepl("Kidney", colnames(df))]
cts_df = df[,c("Geneid", "FetalKidney_438", "FetalKidney_500", "FetalKidney_503", "Kidney_M08", "Kidney_M22", "CortexOfKidney_2181", "CortexOfKidney_6819","FetalKidney_2181","FetalKidney_6819", "Kidney_Daisy")]
myDEA("Kidney", cts_df)

#### wroking for Brain
colnames(df)[grepl("Fetal", colnames(df))]
colnames(df)[grepl("Brain", colnames(df))]
cts_df = df[,c("Geneid", "FetalBrain_410", "FetalBrain_438", "FetalBrain_500", "FetalBrain_503", "BrainCaudalLobe_2181", "BrainCaudalLobe_6819", "Brain_Daisy","Brainstem_2181","Brainstem_6819", "FetalBrain_2181", "FetalBrain_6819", "Cortex_M08", "Cortex_M22")]
myDEA("Brain", cts_df)
