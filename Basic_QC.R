library(data.table)
library(tidyverse)
samples = fread("01_data_preprocess/QC_passed.list", header=F) %>% pull(V1) %>% as.vector

### peaks
mydf = data.frame()
for (s in 1:length(samples)){
    file = paste0("01_data_preprocess/07_peak_called/", samples[s], "_Peaks.bed")
    df = fread(file, header=F) %>% as.data.frame()
    mydf = rbind(mydf, data.frame(Sample=samples[s], NumPeaks=nrow(df)))
    rm(df)
}
mydf %>% separate(., Sample, into=c("Assay", "Tissue", "ID"), sep="_") -> mydf
mydf$Assay = factor(mydf$Assay, levels=c("ATAC", "CTCF", "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3"))
assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999"))
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df[match(colors, color_df$colors),"Assay"])
p=ggplot(mydf) + 
      geom_boxplot(aes(x=Assay, y=NumPeaks/1000, fill=Assay), notch=T, outlier.shape = NA)+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab(expression("# Peaks (×10"^3*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave("Number_of_peaks_by_assays.pdf", p, width=6, height=4.5)
save(mydf, file="Number_of_peaks_by_assays.pdf.Rdata")

p=ggplot(mydf %>% filter(Assay != "ATAC")) + 
      geom_boxplot(aes(x=Assay, y=NumPeaks/1000, fill=Assay), notch=T, outlier.shape = NA)+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab(expression("# Peaks (×10"^3*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave("Number_of_peaks_by_assays.withoutATAC.pdf", p, width=6, height=4.5)


### peaks
mydf = data.frame()
for (s in 1:length(samples)){
    file = paste0("01_data_preprocess/07_peak_called/", samples[s], "_Peaks.bed")
    df = fread(file, header=F) %>% as.data.frame()
    mydf = rbind(mydf, data.frame(Sample=samples[s], NumPeaks=nrow(df)))
    rm(df)
}
mydf$Source = ifelse(grepl("Daisy|2181|6819",mydf$Sample), "Prowse-Wilkins et al.", "FAANG")
mydf %>% separate(., Sample, into=c("Assay", "Tissue", "ID"), sep="_") -> mydf
mydf$Assay = factor(mydf$Assay, levels=c("ATAC", "CTCF", "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3"))
assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999"))
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df[match(colors, color_df$colors),"Assay"])
mydf %>% subset(Assay != "H3K9me3" & Assay != "H3K36me3") -> mydf

p=ggplot(mydf) + 
      geom_boxplot(aes(x=Assay, y=NumPeaks/1000, fill=Source), notch=T, outlier.shape = NA)+
      scale_fill_manual(values=c("#FF4500", "#4169E1"))+
	    xlab("") +
      ylab(expression("# Peaks (×10"^3*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "bottom",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
ggsave("Number_of_peaks_by_assays_comparison.pdf", p, width=5, height=4.5)

library(rstatix)
mydf %>% group_by(Assay) %>% t_test(NumPeaks ~ Source)





samples = fread("01_data_preprocess/QC_passed.list", header=F) %>% pull(V1) %>% as.vector

## spp output format

### spp summary
mydf = data.frame()
for (s in 1:length(samples)){
    file = paste0("01_data_preprocess/06_metrics/", samples[s], ".spp_stats.txt")
    df = fread(file, header=F) %>% as.data.frame()
    names(df)=c("Sample", "numReads", "estFragLen", "corr_estFragLen", "PhantomPeak", "corr_phantomPeak", "argmin_corr", "min_corr", "NSC", "RSC", "QualityTag")
    mydf = rbind(mydf, df)
    rm(df)
}
mydf$Sample = str_replace(mydf$Sample, ".dedup.bam","") %>% str_replace(., "Batch*.","")
mydf %>% separate(., Sample, into=c("Assay", "Tissue", "ID"), sep="_") -> mydf

### color legend
mydf$Assay = factor(mydf$Assay, levels=c("ATAC", "CTCF", "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3"))
assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999"))
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df[match(colors, color_df$colors),"Assay"])

######### plot NSC
mydf[!is.na(mydf$NSC),] -> plot_df
p=ggplot(plot_df) + 
      geom_boxplot(aes(x=Assay, y=NSC, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=NSC), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,5))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab("Normalized strand cross-coefficient\n(NSC)")+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(4,1,1,1), "cm"))
ggsave("all_samples.NSC_byassay.pdf", p, width=6, height=6)
save(mydf, file="all_samples.NSC_byassay.pdf.Rdata")

p=ggplot(plot_df %>% filter(Assay != "ATAC") ) + 
      geom_boxplot(aes(x=Assay, y=NSC, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=NSC), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,5))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab("Normalized strand cross-coefficient\n(NSC)")+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(4,1,1,1), "cm"))
ggsave("all_samples.NSC_byassay.without_ATAC.pdf", p, width=6, height=6)
save(mydf, file="all_samples.NSC_byassay.without_ATAC.pdf.Rdata")


######### plot RSC
mydf[!is.na(mydf$RSC),] -> plot_df
p=ggplot(plot_df) + 
      geom_boxplot(aes(x=Assay, y=RSC, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=RSC), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,5))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab("Relative cross-coefficient\n(RSC)")+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(4,1,1,1), "cm"))
ggsave("all_samples.RSC_byassay.pdf", p, width=6, height=6)
save(mydf, file="all_samples.RSC_byassay.pdf.Rdata")

p=ggplot(plot_df %>% filter(Assay != "ATAC") ) + 
      geom_boxplot(aes(x=Assay, y=RSC, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=RSC), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,5))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab("Relative cross-coefficient\n(RSC)")+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(4,1,1,1), "cm"))
ggsave("all_samples.RSC_byassay.without_ATAC.pdf", p, width=6, height=6)




samples = fread("01_data_preprocess/QC_passed.list", header=F) %>% pull(V1) %>% as.vector
stats=data.frame()
for (s in samples){
  stat_df = fread(paste0("01_data_preprocess/06_metrics/", s,".num_reads.txt"), header=F) %>% as.data.frame
  names(stat_df)="numReads"
  stat_df$Sample = s
  stats=rbind(stats, stat_df)
}
stats$Sample =  str_replace(stats$Sample, "Batch*.", "")
stats %>% separate(., Sample, into=c("Assay", "Tissue", "ID"), sep="_") -> stats
stats$Assay = factor(stats$Assay, levels=c("ATAC", "CTCF", "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3", "Input"))

### color legend
assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF", "Input")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999","#202020"))
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df[match(colors, color_df$colors),"Assay"])

######### plot number of reads
p=ggplot(stats) + 
      geom_boxplot(aes(x=Assay, y=numReads/1000000, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=numReads/1000000), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,100))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab(expression("# usable reads (×10"^6*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(2,1,1,1), "cm"))
ggsave("all_samples.usable_reads.pdf", p, width=6, height=5)
save(stats, file="all_samples.usable_reads.pdf.Rdata")

######### plot number of reads
p=ggplot(stats %>% filter(Assay != "ATAC")) + 
      geom_boxplot(aes(x=Assay, y=numReads/1000000, fill=Assay), notch=F, outlier.shape = NA)+
      geom_jitter(aes(x=Assay, y=numReads/1000000), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,50))+
      scale_fill_manual(values=c(colors))+
	    xlab("") +
      ylab(expression("# usable reads (×10"^6*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "none",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(2,1,1,1), "cm"))
ggsave("all_samples.usable_reads.without_ATAC.pdf", p, width=6, height=5)





samples = fread("01_data_preprocess/QC_passed.list", header=F) %>% pull(V1) %>% as.vector
stats=data.frame()
for (s in samples){
  stat_df = fread(paste0("01_data_preprocess/06_metrics/", s,".num_reads.txt"), header=F) %>% as.data.frame
  names(stat_df)="numReads"
  stat_df$Sample = s
  stats=rbind(stats, stat_df)
}
stats$Sample =  str_replace(stats$Sample, "Batch*.", "")
stats$Source = ifelse(grepl("Daisy|2181|6819",stats$Sample), "Prowse-Wilkins et al.", "FAANG")
stats %>% separate(., Sample, into=c("Assay", "Tissue", "ID"), sep="_") -> stats
stats$Assay = factor(stats$Assay, levels=c("ATAC", "CTCF", "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3", "Input"))

### color legend
assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF", "Input")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999","#202020"))
colors = as.vector(color_df$colors)
names(colors)=as.vector(color_df[match(colors, color_df$colors),"Assay"])

stats %>% subset(Assay != "H3K9me3" & Assay != "H3K36me3") -> stats
######### plot number of reads
p=ggplot(stats) + 
      geom_boxplot(aes(x=Assay, y=numReads/1000000, fill=Source), notch=F, outlier.shape = NA)+
      #geom_jitter(aes(x=Assay, y=numReads/1000000), shape=16, size=1.5, position=position_jitter(0.2))+
      ylim(c(0,100))+
      scale_fill_manual(values=c("#FF4500", "#4169E1"))+
	    xlab("") +
      ylab(expression("# usable reads (×10"^6*")"))+
      theme_classic(base_size = 20, base_line_size = 1)+
      theme(legend.position = "bottom",
            axis.text.x = element_text(color="black", size=15, angle = 45, hjust=1),
            axis.text.y = element_text(color="black", size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(2,1,1,1), "cm"))
ggsave("Usable_reads_comparison.pdf", p, width=5, height=5)
