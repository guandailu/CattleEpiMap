library(GenomicRanges)
library(data.table)
library(tidyverse)
library(rGREAT)
ARGS <- commandArgs(trailingOnly = TRUE)
state = ARGS[1]
# tissue_specific_RE_analysis.R
out_dir=paste0(state, "/Tisssue_specific_rGREAT")
if (!dir.exists(out_dir)){
    dir.create(out_dir)
}
files=list.files(paste0(state, "/consensus_state/"),".bed", full.names=F)
for (f in 1:length(files)){
  df = fread(paste0(state, "/consensus_state/", files[f]), header=F) %>% as.data.frame
  names(df)=c("chr", "start", "end", "isREs")
  df$Tissue = str_replace(files[f], ".bed","")
  if (f == 1){
    mydf = df
  }else{
    mydf = rbind(mydf, df)
  }
}
mydf %>% group_by(Tissue) %>% summarize(num_all=sum(isREs)) -> plot_df1

## define colors
meta = read.table("01_data_preprocess/sample_clustering/samples_meta.txt", header=T, comment.char="@", sep="\t", stringsAsFactors=T) %>% as.data.frame
colors = paste0("#",meta$Tissue_colors)
names(colors)=as.vector(meta$Tissue)
colors = colors[!duplicated(colors)]
## plot all RE per tissue
p = ggplot(plot_df1)+
    geom_bar(aes(x=reorder(Tissue, num_all), y = num_all, fill=Tissue),stat="identity")+
    coord_flip()+
    xlab("") +
    ylab("#  elements")+
    scale_fill_manual(values=c(colors))+
    theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0(state, "/Num_", state,"_by_tissue.pdf"), p, width=4, height=6)

files=list.files(paste0(state, "/consensus_state/"),paste0("_specific_",state,".txt"), full.names=F)
for (f in 1:length(files)){
  df = fread(paste0(state, "/consensus_state/", files[f]), header=F) %>% as.data.frame
  names(df)=c("chr", "start", "end", "isREs")
  df$Tissue = str_replace(files[f], paste0("_specific_",state,".txt"),"")
  if (f == 1){
    mydf = df
  }else{
    mydf = rbind(mydf, df)
  }
}
mydf %>% group_by(Tissue) %>% summarize(num_spefic=sum(isREs)) -> plot_df2
## plot tissue specific RE per tissue
p = ggplot(plot_df2)+
    geom_bar(aes(x=reorder(Tissue, num_spefic), y = num_spefic, fill=Tissue),stat="identity")+
    coord_flip()+
    xlab("") +
    ylab("# tissue specific EnhA")+
    scale_fill_manual(values=c(colors))+
    theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(paste0(state, "/Num_tissue_specific_", state,".pdf"), p, width=4, height=6)


## plot both tissue specific RE and all RE per tissue
plot_df = merge(plot_df1, plot_df2, by = "Tissue")
ylim.prim <- c(0, max(plot_df$num_all))
ylim.sec <- c(0, max(plot_df$num_spefic))
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
p = ggplot(plot_df)+
    geom_bar(aes(x=reorder(Tissue, num_all), y = num_all, fill=Tissue),stat="identity")+
    geom_point(aes(x=Tissue, y= a + num_spefic*b), color="black",shape=19, size=3)+
    geom_line(aes(x=Tissue, y = a + num_spefic*b, group=1), color="black", linewidth=0.8)+
    scale_y_continuous(name = "# regulatory elements",sec.axis = sec_axis(trans=~(.-a)/b, name="# tissue specific"))+
    coord_flip()+
    xlab("") +
    scale_fill_manual(values=c(colors))+
    theme_classic(base_size = 15, base_line_size = 0.8)+
    theme(legend.position="none", axis.text.x = element_text(color="black", size=10),
      axis.text.y.left = element_text(color="black", size=10),
      axis.title.y.right= element_text(color="red", size=10),
      axis.text.y.right= element_text(color="red", size=10),
      axis.ticks.y.right= element_line(color="red"),
      axis.line.y.right= element_line(color="red"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
ggsave(paste0(state, "/Num_tissue_all_and_specific_", state,".pdf"), p, width=4, height=6)

# rGreat enrichment analysis
files=list.files(paste0(state, "/consensus_state/"), paste0("_specific_",state,".txt"), full.names=F)
for (f in 1:length(files)){
  df = fread(paste0(state, "/consensus_state/", files[f]), header=F) %>% as.data.frame
  names(df)=c("chr", "start", "end", "isREs")
  df$Tissue = str_replace(files[f], paste0("_specific_",state,".txt"),"")
  cat("****number of ", files[f]," specific REs: ", nrow(df), "\n")
  gr=GRanges(seqnames=c(df$chr),
             ranges=IRanges(start=c(as.integer(df$start)),
                            end=c(as.integer(df$end))),
             strand=c(rep("*", nrow(df)))
  )
  job=great(gr, "GO:BP", "bosTau9")
  prefix=str_replace(files[f], paste0("_specific_",state,".txt"),"")
  pdf(paste0(out_dir, "/", prefix,".GO_BP_plotVolcano.pdf"), width=8, height=6)
  plotVolcano(job)
  dev.off()

  pdf(paste0(out_dir, "/", prefix,".GO_BP_plotRegionGeneAssociations.pdf"), width=12, height=6)
  plotRegionGeneAssociations(job)
  dev.off()

  res_df=getRegionGeneAssociations(job) %>% as.data.frame
  jobres = job@table %>% arrange(p_value)
  apply(res_df,2,as.character) -> res_df
  write.table(res_df, paste0(out_dir, "/", prefix, ".rGreat_dist2gene.txt"), row.names=F, sep="\t", quote=F)
  write.table(jobres, paste0(out_dir,"/", prefix, ".rGreat_enrichGO_BP.txt"), row.names=F, sep="\t", quote=F)
}




library(data.table)
library(tidyverse)

files=list.files(path ="E4/Tisssue_specific_HOMER", pattern="knownResults.txt$", recursive = TRUE, full.names = TRUE)

mydf = data.frame()
for (i in 1:length(files)){
df = fread(files[i]) %>% as.data.frame
names(df)[1:5]=c("MotifName", "Consensus", "pval", "logp", "qval")
df=df[!grepl("SeqBias", df$MotifName),]
subdf = subset(df, qval < 0.05)
if (nrow(subdf) > 0){
  subdf = subdf 
}else{
  #subdf = subset(df, pval < 0.05)
  #if (nrow(subdf) > 0){
  #  subdf = subdf
  #}else{
  #  subdf = df[1,]
  #}
  subdf = data.frame()
}
if (nrow(subdf) > 0){
subdf = subdf %>% mutate(logp=-log10(qval)) %>% select(MotifName, logp) %>% mutate(Tissue=str_replace(files[i], "/knownResults.txt", "") %>% str_replace(., "E4/Tisssue_specific_HOMER/", ""))
mydf = rbind(mydf, subdf)
}
}
mydf[sapply(mydf, is.infinite)] <- NA
mydf[is.na(mydf)] <- max(mydf$logp, na.rm=T)

#mydf$logp[mydf$logp > 30] <- 30
mydf = mydf %>% group_by(MotifName, Tissue) %>% summarize(logp = max(logp, na.rm = TRUE)) %>% spread(Tissue, logp, fill=0) %>% as.data.frame
rownames(mydf)=mydf$MotifName
mydf$MotifName = NULL


meta = read.table("/group/zhougrp4/dguan/BovineFAANG/01_data_preprocess/sample_clustering/samples_meta.txt", header=T, comment.char="@", sep="\t", stringsAsFactors=T) %>% as.data.frame
colors = paste0("#",meta$Tissue_colors)
names(colors)=as.vector(meta$Tissue)
colors = colors[!duplicated(colors)]
ann_colors=list(Tissue=colors[names(colors) %in% names(mydf)])
annotation_row = data.frame(Tissue=unique(names(colors)))
rownames(annotation_row)=annotation_row$Tissue
annotation_row=subset(annotation_row, Tissue %in% names(mydf))
colnames(mydf)[colnames(mydf) == "LymphNodes"] <- "LymphNode"
library(pheatmap)
pdf("HOMER_enrichment_all.pdf", width=6, height=8)
pheatmap(as.matrix(mydf),
         color = colorRampPalette(c("#FFFFFF", "#FFC000", "#C73653"))(50),
         annotation_col = annotation_row,
         annotation_colors=ann_colors,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = T,
         border_color = NA,
         #display_numbers = matrix(ifelse(mydf > 1, "*", ""), nrow(mydf)),
         #angle_col = 90
        )
dev.off()


pdf("HOMER_enrichment_all.long.pdf", width=24, height=120)
pheatmap(as.matrix(mydf),
         color = colorRampPalette(c("#FFFFFF", "#FFC000", "#C73653"))(50),
         annotation_col = annotation_row,
         annotation_colors=ann_colors,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         #display_numbers = matrix(ifelse(mydf > 1, "*", ""), nrow(mydf)),
         #angle_col = 90
        )
dev.off()




########## top 5
mydf = data.frame()
for (i in 1:length(files)){
df = fread(files[i]) %>% as.data.frame
names(df)[1:5]=c("MotifName", "Consensus", "pval", "logp", "qval")
df=df[!grepl("SeqBias", df$MotifName),]
subdf = df[1:5,]
if (nrow(subdf) > 0){
  subdf = subdf 
}else{
  #subdf = subset(df, pval < 0.05)
  #if (nrow(subdf) > 0){
  #  subdf = subdf
  #}else{
  #  subdf = df[1,]
  #}
  subdf = data.frame()
}
if (nrow(subdf) > 0){
subdf = subdf %>% mutate(logp=-log10(pval)) %>% select(MotifName, logp) %>% mutate(Tissue=str_replace(files[i], "/knownResults.txt", "") %>% str_replace(., "E4/Tisssue_specific_HOMER/", ""))
mydf = rbind(mydf, subdf)
}
}
mydf[sapply(mydf, is.infinite)] <- NA
mydf[is.na(mydf)] <- max(mydf$logp, na.rm=T)

mydf$logp[mydf$logp > 20] <- 20
mydf = mydf %>% group_by(MotifName, Tissue) %>% summarize(logp = max(logp, na.rm = TRUE)) %>% spread(Tissue, logp, fill=0) %>% as.data.frame
rownames(mydf)=mydf$MotifName
mydf$MotifName = NULL


meta = read.table("/group/zhougrp4/dguan/BovineFAANG/01_data_preprocess/sample_clustering/samples_meta.txt", header=T, comment.char="@", sep="\t", stringsAsFactors=T) %>% as.data.frame
colors = paste0("#",meta$Tissue_colors)
names(colors)=as.vector(meta$Tissue)
colors = colors[!duplicated(colors)]
ann_colors=list(Tissue=colors[names(colors) %in% names(mydf)])
annotation_row = data.frame(Tissue=unique(names(colors)))
rownames(annotation_row)=annotation_row$Tissue
annotation_row=subset(annotation_row, Tissue %in% names(mydf))

library(pheatmap)
pdf("HOMER_enrichment_all.top5.pdf", width=6, height=8)
pheatmap(as.matrix(mydf),
         color = colorRampPalette(c("#FFFFFF", "#FFC000", "#C73653"))(50),
         annotation_col = annotation_row,
         annotation_colors=ann_colors,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = T,
         border_color = NA,
         #display_numbers = matrix(ifelse(mydf > 1, "*", ""), nrow(mydf)),
         #angle_col = 90
        )
dev.off()


pdf("HOMER_enrichment_all.long.top5.pdf", width=24, height=120)
pheatmap(as.matrix(mydf),
         color = colorRampPalette(c("#FFFFFF", "#FFC000", "#C73653"))(50),
         annotation_col = annotation_row,
         annotation_colors=ann_colors,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         #display_numbers = matrix(ifelse(mydf > 1, "*", ""), nrow(mydf)),
         #angle_col = 90
        )
dev.off()