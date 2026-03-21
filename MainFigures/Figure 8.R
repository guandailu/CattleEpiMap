library(data.table)
library(tidyverse)

df = fread("03_human2cattle_CREs/CRE.sum.txt", header=F)
names(df)=c("CRE", "Type", "num")

p=ggplot(df)+
    geom_bar(aes(x=factor(CRE, levels=c(paste0("E", 14:1))), y =num, fill=factor(Type, levels=c("sfCRE", "sdCRE", "soCRE", "ssCRE"))), position="fill", stat="identity")+
    coord_flip()+
    ylab("Proportion of CREs")+
    xlab("")+
    scale_fill_manual(values=c(colors))+
    theme_classic(base_size = 20, base_line_size = 0.8)+
    theme(legend.position="bottom",legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=15, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
          guides(fill=guide_legend(nrow=1))
ggsave(paste0("03_human2cattle_CREs/CRE.sum.pdf"), p , width=3, height=6)


library(data.table)
library(tidyverse)

df_list=list()
for (type in c("sfCRE", "sdCRE", "soCRE", "ssCRE")){
  df = fread(paste0("06_CRE_toTSS/E5.", type,".toTSS.out"), header=F)
  df$Type=type
  df_list[[type]]=df
}
rbindlist(df_list) -> df

df$bins = cut(df$V12, breaks=c(min(df$V12), -1000000, -500000, -100000, -50000, 0, 50000, 100000, 500000,  max(df$V12)), labels=c(-1000, -500, -100, -50, 0, 50, 100, 500, 1000), include.lowest= T)

df %>% group_by(Type, bins) %>% summarize(num=length(unique(paste0(V1, "_",V2,"_",V3)))) -> plot_df
plot_df %>% group_by(Type) %>% mutate(perc=num / sum(num)) -> plot_df
library("RColorBrewer")

p=ggplot(plot_df, aes(x=bins, y=perc, fill=factor(Type, levels=c("sfCRE", "sdCRE", "soCRE", "ssCRE"))))+
    geom_bar(stat="identity", position=position_dodge())+
    #geom_point()+
    #geom_line()+
    ylab("Proportion of EnhA")+
    xlab("")+
    scale_fill_manual(values=c(colors))+
    theme_classic(base_size = 20, base_line_size = 0.8)+
    theme(legend.position="bottom",legend.title = element_blank(),
          axis.text.x = element_text(color="black", size=15, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color="black", size=15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
          guides(fill=guide_legend(nrow=1))
ggsave(paste0("06_CRE_toTSS/EnhA_toTSS_by_conservation_types.pdf"), p , width=4, height=4)


