library(data.table)
library(tidyverse)
bed_files=list.files("find_markers", pattern="Markers..*.bed$")
for (b in 1:length(bed_files)){
  df = fread(paste0("find_markers/", bed_files[b]), header=T) %>% as.data.frame %>% filter(direction == "U")
  if (nrow(df) > 50){
    df = df[1:50,]
  }
  if (b == 1){
    mydf = df
  }else{
    mydf = rbind(mydf, df)
  }
}
tiss_spec_cpgs=paste(mydf$`#chr`, mydf$start, mydf$end, mydf$startCpG, mydf$endCpG, sep="_") %>% unique
tiss_spec_samples=fread("groups_tissues_DMR.csv", header=T) %>% as.data.frame %>% mutate(Sample=str_replace(name, ".sorted","")) %>% pull(Sample)
counts=fread("beta_to_table/BovineFAANG_WGBS_segments.counts.txt", header=T) %>% as.data.frame
rownames(counts)=paste(counts$chr, counts$start, counts$end, counts$startCpG, counts$endCpG, sep="_")
counts$chr=counts$start=counts$end=counts$startCpG=counts$endCpG=NULL
names(counts)=str_replace(colnames(counts), ".sorted","")
dim(counts)
counts[tiss_spec_cpgs,tiss_spec_samples] -> counts
counts[is.na(counts)] <- 0

library(pheatmap)
library(RColorBrewer)

pdf("find_markers/tissue_specific_methy_fragments.2.pdf", width=6, height=7)
pheatmap(as.matrix(counts),
         color = colorRampPalette(c("#003366", "white", "#FFFFCC"))(100),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = T,
       fontsize =8
        )
dev.off()






library(data.table)
library(tidyverse)
### data load

m="ATAC"

files=list.files(path="regulatory_map/ATAC", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/ATAC/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 100, 200))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/ATAC/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)



library(data.table)
library(tidyverse)
### data load

m="H3K27ac"

files=list.files(path="regulatory_map/H3K27ac/", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/H3K27ac/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
head(mydf)
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 75, 150))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/H3K27ac/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)




library(data.table)
library(tidyverse)
### data load

m="H3K4me1"

files=list.files(path="regulatory_map/H3K4me1/", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/H3K4me1/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
head(mydf)
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 75, 150))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/H3K4me1/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)






library(data.table)
library(tidyverse)
### data load

m="H3K4me3"

files=list.files(path="regulatory_map/H3K4me3/", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/H3K4me3/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
head(mydf)
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 75, 150))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/H3K4me3/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)




library(data.table)
library(tidyverse)
### data load

m="H3K27me3"

files=list.files(path="regulatory_map/H3K27me3/", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/H3K27me3/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
head(mydf)
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 75, 150))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/H3K27me3/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)




library(data.table)
library(tidyverse)
### data load

m="CTCF"

files=list.files(path="regulatory_map/CTCF/", pattern="Cerebellum_hypo.*.csv$")

for (f in 1:length(files)){
  file=paste0("regulatory_map/CTCF/", files[f])
  if (file.size(file) != 0){
    df = fread(file, header=T, skip=1) %>% as.data.frame
    df %>% gather(., Position, Enrichment, 3:ncol(df)) -> df
    names(df)[1:2]=c("Sample", "File")
    df$Sample = str_replace(files[f],".matrix.mat.csv","") %>% str_replace(., ".*_in_","")
    df$Mark = str_replace(files[f],"_signals_in_.*csv$","") %>% str_replace(., ".*_hypo\\.","")
    df$Enrichment = as.numeric(as.character(df$Enrichment))
    df$Position = as.integer(as.character(df$Position))
    if (f ==1){
      mydf = df
    }else{
      mydf = rbind(mydf, df)
    }
  }
}
na.omit(mydf) -> mydf
separate(mydf, Sample, into=c("Tissue", "ID"), sep="_", remove=F) -> mydf
head(mydf)
mydf %>% group_by(Mark, Tissue, Position) %>% summarize(tissEnrichment = mean(Enrichment)) -> mydf
mydf %>% mutate(Tissue = case_when(Tissue == "Cerebellum" ~ "Cerebellum", .default = "Background")) -> plotdf
plotdf %>% group_by(Mark, Tissue, Position) %>% summarize(Enrichment = mean(tissEnrichment)) -> plotdf

scale_min_1 <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  scaled <- ((x - min_val) / (max_val - min_val)) * (max_val - 1) + 1
  return(scaled)
}
plotdf %>% group_by(Tissue) %>% mutate(Enrichment = scale_min_1(Enrichment)) -> plotdf

p=ggplot(data=plotdf, aes(x=Position, y=Enrichment, color=Tissue)) +
    geom_line()+
    #ylim(0, 80)+
    #facet_wrap(~factor(Mark), nrow =1)+
    scale_x_continuous("", c(1, 75, 150))+
    scale_color_manual(values= c("#E0E0E0", "#4C9900"))+
    #scale_color_manual(values= c(colors))+
    xlab("") +
    ylab("Enrichment")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position.inside=c(0.2, 0.8), 
          legend.title=element_text(), 
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0("regulatory_map/CTCF/Cerebellum.hypo_regions_regulatory_map_", m,".pdf"), p, width=6, height=4)



#### histone modification
# computeMatrix reference-point --referencePoint center --binSize 1000 -S ${bw_file} -R ${tmp_file} -o regulatory_map/${mark}/${tissue}_hypo.${mark}_signals_in_${sample}.matrix.mat.gz -b 75000 -a 75000 --numberOfProcessors 12 --skipZeros --scale 1


# ATAC
# computeMatrix reference-point --referencePoint center --binSize 100 -S ${bw_file} -R ${tmp_file} -o regulatory_map/${mark}/${tissue}_hypo.${mark}_signals_in_${sample}.matrix.mat.gz -b 10000 -a 10000 --numberOfProcessors 12 --skipZeros --scale 2
