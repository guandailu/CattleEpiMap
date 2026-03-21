library(rjson)
library(tidyverse)

files=system("ls ATAC_bams/*/chrombpnet_model/evaluation/chrombpnet_metrics.json", intern=T)
mydf=data.frame()
for (f in 1:length(files)){
  tissue=str_replace(files[f], "ATAC_bams/", "") %>% str_replace(., "/chrombpnet_model/evaluation/chrombpnet_metrics.json","")
  df = fromJSON(file =files[f]) %>% as.data.frame %>% mutate(Tissue=tissue)
  mydf = rbind(mydf, df)
}
write.table(mydf, "chrombpnet_metrics.summary.txt", row.names=F, sep="\t", quote=F)
mydf[mydf$Tissue=="Myoblast","Tissue"] = "Myoblasts"

## define colors
meta = read.table("01_data_preprocess/sample_clustering/samples_meta.txt", header=T, comment.char="@", sep="\t", stringsAsFactors=T) %>% as.data.frame
colors = paste0("#",meta$Tissue_colors)
names(colors)=as.vector(meta$Tissue)
colors = colors[!duplicated(colors)]
colors = colors[mydf$Tissue]

p=ggplot(data=mydf, aes(x=reorder(Tissue, -counts_metrics.peaks.spearmanr), y=counts_metrics.peaks.spearmanr, fill=Tissue)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c(colors))+
    xlab("") +
    ylab("Spearman's rho")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_counts_metrics.spearmanr.pdf", p, width=7, height=4.2)


p=ggplot(data=mydf, aes(x=reorder(Tissue, -counts_metrics.peaks.pearsonr), y=counts_metrics.peaks.pearsonr, fill=Tissue)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c(colors))+
    xlab("") +
    ylab("Pearson's r")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_counts_metrics.pearsonr.2.pdf", p, width=7, height=4.2)


p=ggplot(data=mydf %>% select(counts_metrics.peaks.pearsonr, Tissue) %>% mutate(type="pearsonr"), aes(x=type,y=counts_metrics.peaks.pearsonr)) +
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=Tissue),shape=16, position=position_jitter(0.2), size=3)+
    scale_y_continuous(name=c(seq(0.3, 0.9, 0.1)), breaks=c(seq(0.3, 0.9, 0.1)), limits=c(0.3, 0.9))+
    #geom_point()+
    scale_color_manual(values=c(colors))+
    coord_flip()+
    xlab("") +
    ylab("Pearson's r")+
    theme_linedraw(base_size = 20, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=15, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_counts_metrics.pearsonr.boxplot.pdf", p, width=7, height=2)



p=ggplot(data=mydf %>% select(counts_metrics.peaks.pearsonr, Tissue) %>% mutate(type="pearsonr"), aes(x=type,y=counts_metrics.peaks.pearsonr)) +
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=Tissue),shape=16, position=position_jitter(0.2), size=3)+
    scale_y_continuous(name=c(seq(0.3, 0.9, 0.1)), breaks=c(seq(0.3, 0.9, 0.1)), limits=c(0.3, 0.9))+
    #geom_point()+
    scale_color_manual(values=c(colors))+
    coord_flip()+
    xlab("") +
    ylab("Pearson's r")+
    theme_classic(base_size = 20, base_line_size = 0.5)+
    theme(
          axis.text.x = element_text(color="black", size=15, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_counts_metrics.pearsonr.boxplot.with_lengend.pdf", p, width=10, height=8)



p=ggplot(data=mydf, aes(x=reorder(Tissue, -counts_metrics.peaks.mse), y=counts_metrics.peaks.mse, fill=Tissue)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c(colors))+
    xlab("") +
    ylab("MSE")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_counts_metrics.mse.pdf", p, width=7, height=4.2)



p=ggplot(data=mydf, aes(x=reorder(Tissue, -profile_metrics.peaks.median_jsd), y=profile_metrics.peaks.median_jsd, fill=Tissue)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c(colors))+
    xlab("") +
    ylab("median_jsd")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_profile_metrics.median_jsd.pdf", p, width=7, height=4.2)


p=ggplot(data=mydf, aes(x=reorder(Tissue, -profile_metrics.peaks.median_norm_jsd), y=profile_metrics.peaks.median_norm_jsd, fill=Tissue)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=c(colors))+
    xlab("") +
    ylab("median_norm_jsd")+
    theme_classic(base_size = 15, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=10, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_profile_metrics.median_norm_jsd.pdf", p, width=7, height=4.2)


mydf %>% gather(metrics, value, counts_metrics.peaks.spearmanr:profile_metrics.peaks.median_norm_jsd) %>% group_by(metrics)  %>% summarize(val=median(value)) -> df
library(RColorBrewer)
p=ggplot(data=df, aes(x=metrics, y=val, fill=metrics)) +
    geom_bar(stat="identity")+
    scale_fill_brewer(palette="Paired")+
    xlab("") +
    ylab("Evaluation metrics")+
    theme_classic(base_size = 20, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=15, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_evaluation_metrics.pdf", p, width=4, height=7)



mydf %>% gather(metrics, value, counts_metrics.peaks.spearmanr:profile_metrics.peaks.median_norm_jsd) -> df
library(RColorBrewer)
p=ggplot(data=df, aes(x=metrics, y=value)) +
    geom_boxplot()+
    geom_jitter(aes(color=Tissue), shape=16, position=position_jitter(0.2))+
    scale_color_manual(values=c(colors))+
    xlab("") +
    ylab("Evaluation metrics")+
    theme_classic(base_size = 20, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=15, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_evaluation_metrics.boxplot.pdf", p, width=4, height=7)



df %>% subset(metrics != "counts_metrics.peaks.mse") -> df2
p=ggplot(data=df2, aes(x=metrics, y=value)) +
    geom_boxplot()+
    geom_jitter(aes(color=Tissue), shape=16, position=position_jitter(0.2))+
    scale_color_manual(values=c(colors))+
    xlab("") +
    ylab("Evaluation metrics")+
    theme_classic(base_size = 20, base_line_size = 0.5)+
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=15, angle=90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chrombpnet_evaluation_metrics.boxplot.subset.pdf", p, width=4, height=7)




> df %>% subset(metrics == "counts_metrics.peaks.pearsonr") %>% mutate(value = round(value, 2)) %>% summarise(as_tibble_row(quantile(value)))
    0%  25%  50% 75% 100%
1 0.35 0.72 0.76 0.8 0.84



subset(df, metrics == "counts_metrics.peaks.pearsonr") %>% mutate(value = round(value, 2)) %>% filter(value >= 0.72) %>% arrange(value) %>% select(Tissue, value) %>% rename(pearsonr=value) ->out_df
write.table(out_df, "final_tissues_with_chrombpnet_models.txt", sep="\t", col.names=F, quote=F, row.names=F)



library(rhdf5)
library(data.table)
library(tidyverse)
library(ggseqlogo)
library(pheatmap)

tissues <- fread("final_tissues_with_chrombpnet_models.txt", header=F) %>% pull(V1)

transform_to_pwm <- function(mat){
    mat[mat < 0] <- 0
    mat <- mat + 1e-6
    pwm <- t(apply(mat, 1, function(x) x / sum(x)))
    return(pwm)
}


all_motifs <- list()
for (ts in tissues){
    h5file <- paste0("ATAC_bams/", ts,"/chrombpnet_motifs/", ts,".modisco_results.h5")
    motif_paths <- h5ls(h5file, recursive=TRUE) %>% 
                   filter(grepl("contrib_scores", name)) %>% 
                   filter(!grepl("subpattern|seqlets", group)) %>%
                   pull(group) %>%
                   paste0("/contrib_scores")

    motifs_list <- lapply(motif_paths, function(path) {
        h5read(h5file, path)
    })

    motifs_pwm <- lapply(motifs_list, transform_to_pwm)
    motif_names <- paste0(ts, "_pattern_", seq_along(motifs_pwm)-1)
    names(motifs_pwm) <- motif_names

    all_motifs <- c(all_motifs, motifs_pwm)
}

max_len <- max(sapply(all_motifs, nrow))
pad_pwm <- function(pwm, target_len){
    if(nrow(pwm) < target_len){
        pad_len <- target_len - nrow(pwm)
        pad <- matrix(rep(0.25, pad_len*4), ncol=4)
        pwm <- rbind(pwm, pad)
    }
    return(pwm)
}
all_motifs_padded <- lapply(all_motifs, pad_pwm, target_len=max_len)
motif_mat <- do.call(rbind, lapply(all_motifs_padded, as.vector))
rownames(motif_mat) <- names(all_motifs_padded)
similarity_matrix <- cor(t(motif_mat))
    
system("mkdir -p merged_motifs")
pdf("merged_motifs/merged_motifs.similarity_matrix.pdf", width=8, height=12)                     
pheatmap(similarity_matrix, fontsize=7, main="Motif Similarity")
dev.off()
                   
distance_matrix <- as.dist(1 - similarity_matrix)
hc <- hclust(distance_matrix, method="average")

pdf("merged_motifs/merged_motifs.clustering.pdf", width=24, height=8)                     
plot(hc, cex=0.5, main="Hierarchical Clustering of Motifs")
dev.off()


clusters <- cutree(hc, h=0.1)
cluster_df <- data.frame(
    motif = names(clusters),
    cluster = clusters
)

# Merge motifs within each cluster
merged_motifs <- list()
motif_mapping <- data.frame()

for(cluster_id in unique(clusters)){
    motif_indices <- which(clusters == cluster_id)
    motif_names <- names(all_motifs_padded)[motif_indices]
    motif_matrices <- all_motifs_padded[motif_indices]

    merged_pwm <- Reduce("+", motif_matrices) / length(motif_matrices)
    merged_name <- paste0("merged_cluster_", cluster_id)
    merged_motifs[[merged_name]] <- merged_pwm

    motif_mapping <- rbind(motif_mapping,
                           data.frame(merged_motif = merged_name,
                                      original_motif = motif_names))
}

# Visualize merged motifs
for(motif_name in names(merged_motifs)){
    pwm <- as.data.frame(merged_motifs[[motif_name]])
    rownames(pwm)=c("A", "C", "G", "T")
    p <- ggseqlogo(as.matrix(pwm), method='custom') +
        ggtitle(motif_name) + theme_minimal() +
        theme(axis.text.x = element_blank(),
          axis.text.y = element_text(color="black", size=5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    ggsave(paste0("merged_motifs/",motif_name, ".pdf"), p, width=3, height=1.3)
}
fwrite(motif_mapping, "merged_motifs/merged_motifs_mapping_across_tissues.txt", row.names=F, sep="\t", quote=F)        
                   
                   
# Save merged motifs in MEME format
write_meme <- function(pwm_list, filename = "merged_motifs.meme") {
    fileConn <- file(filename, "w")
    writeLines("MEME version 4\n", fileConn)
    writeLines("ALPHABET= ACGT\n", fileConn)
    writeLines("strands: + -\n", fileConn)
    writeLines("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n", fileConn)
    for (i in seq_along(pwm_list)) {
        pwm <- t(pwm_list[[i]])
        motif_name <- names(pwm_list)[i]
        motif_width <- nrow(pwm)
        writeLines(sprintf("MOTIF %s", motif_name), fileConn)
        writeLines(sprintf("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0", motif_width), fileConn)
        for (j in 1:motif_width) {
            line <- sprintf("%.4f %.4f %.4f %.4f", pwm[j,1], pwm[j,2], pwm[j,3], pwm[j,4])
            writeLines(line, fileConn)
        }
        writeLines("", fileConn) 
    }
    close(fileConn)
}

write_meme(merged_motifs, "merged_motifs/final_merged_motifs.meme")
