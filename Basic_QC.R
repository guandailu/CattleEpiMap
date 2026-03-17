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

