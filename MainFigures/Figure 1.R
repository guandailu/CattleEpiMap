library(data.table)
library(tidyverse)
library(ComplexHeatmap)

df = fread("all_samples.txt", header=T) %>% as.data.frame
df$SampleID = str_replace(df$SampleID, "Batch.*$","")
df$EID = str_replace(df$EID, "Batch.*$","")
df$ID = str_replace(df$ID, "Batch.*$","")

df %>% mutate(value=case_when(DataType=="Experimental" ~ "2", DataType=="Imputed" ~"1", .default="0")) %>% select(EID, Mark, value) %>% spread(EID, value, fill="0") -> plot_df
rownames(plot_df)=plot_df$Mark
plot_df$Mark=NULL
plot_df = plot_df[c("RNASeq", "H3K9me3", "H3K36me3", "WGBS", "ATAC", "CTCF", "H3K27me3", "H3K27ac", "H3K4me3", "H3K4me1"),]
apply(plot_df, 2, as.integer) %>% colSums %>% sort(decreasing=T) -> sorted_colnames
names(sorted_colnames) -> sorted_colnames
plot_df = plot_df[, c(sorted_colnames)]


assays=c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "CTCF", "WGBS", "RNASeq")
color_df=data.frame(Assay=assays, colors=c("#53B7E6", "#FAE54C", "#EA3323", "#F4BB40", "#8A999A", "#928E51", "#53B651", "#E39999", "#9370DB", "#4682B4"))
rownames(color_df)=color_df$Assay
color_df = as.data.frame(color_df[rownames(plot_df),])

mark_df = df %>% group_by(Mark) %>% summarize(num=n()) %>% as.data.frame
rownames(mark_df)=mark_df$Mark
mark_df=as.data.frame(mark_df[rownames(plot_df),])

##### color row
mark_colors = as.vector(color_df$colors)
names(mark_colors)=as.vector(color_df$Assay)
row_ha = HeatmapAnnotation(mark=mark_df$Mark,
                           num=anno_barplot(mark_df$num), 
                           col = list(mark = mark_colors),
                           which = "row"
                       )
##### color column
tissue_color_df = df[,c("Tissue", "Tissue_colors")] %>% distinct
tissue_colors=as.vector(paste0("#",tissue_color_df$Tissue_colors))
names(tissue_colors)=as.vector(tissue_color_df$Tissue)

development_color_df = df[,c("Development", "Development_colors")] %>% distinct
development_colors = as.vector(paste0("#",development_color_df$Development_colors))
names(development_colors)=as.vector(development_color_df$Development)

ref_color_df = df[,c("Reference", "Reference_colors")] %>% distinct
ref_colors = as.vector(paste0("#",ref_color_df$Reference_colors))
names(ref_colors)=as.vector(ref_color_df$Reference)

state_sample_df = fread("/group/zhougrp4/dguan/BovineFAANG/03_chromatin_ann/FourMarks/ChromatinStates/sample.txt", header=F)
df$CoreAnn = ifelse(df$EID %in% state_sample_df$V1, "CoreAnn","Not")
CoreAnn_colors=c("#6B8E23", "#FFFFFF")
names(CoreAnn_colors)=c("CoreAnn", "Not")

state_sample_df = fread("/group/zhougrp4/dguan/BovineFAANG/03_chromatin_ann/SixMarks/ChromatinStates/samples.txt", header=F)
df$ExpandAnn = ifelse(df$EID %in% state_sample_df$V1, "ExpandAnn","Not")
ExpandAnn_colors=c("#FF4500", "#FFFFFF")
names(ExpandAnn_colors)=c("ExpandAnn", "Not")

col_ha = HeatmapAnnotation(Tissue=df[match(colnames(plot_df), df$EID),"Tissue"],
                           Development=df[match(colnames(plot_df), df$EID),"Development"],
                           Source=df[match(colnames(plot_df), df$EID),"Reference"],
                           CoreAnn=df[match(colnames(plot_df), df$EID),"CoreAnn"],
                           ExpandAnn=df[match(colnames(plot_df), df$EID),"ExpandAnn"],
                           col = list(Tissue = tissue_colors, Development = development_colors, Source=ref_colors, CoreAnn=CoreAnn_colors, ExpandAnn=ExpandAnn_colors),
                           which = "column"
                       )




colors = c("#FF7F50", "#B22222", "#BFBFBF")
names(colors) = c("1", "2", "0")
ht <- Heatmap(as.matrix(plot_df), 
              name = " ", 
              col = colors, 
              right_annotation = row_ha, 
              bottom_annotation = col_ha, 
              column_title = "", 
              row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 0), 
              row_names_side = c("left"))
pdf(file="all_samples.txt.pdf", width=16, height=3)
draw(ht)
dev.off()


epigenome_df = data.frame(epigenome=colnames(plot_df)) %>% separate(epigenome, into=c("Tissue", "ID"), remove=F, sep="_")
fwrite(epigenome_df, "all_epigneomes.txt", sep="\t", quote=F)