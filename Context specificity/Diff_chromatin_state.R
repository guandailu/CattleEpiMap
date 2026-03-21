library(VennDiagram)
library(tidyverse)
library(data.table)
library(RColorBrewer)

ARGS <- commandArgs(trailingOnly = TRUE)
tissue = ARGS[1]
state = ARGS[2]

fetalBed=list.files("DiffStates/", pattern=paste0("^Fetal", tissue, "(.)*.", state, ".bed"))
for (f in 1:length(fetalBed)){
    fetaldf = fread(paste0("DiffStates/", fetalBed[f]), header=F) %>% as.data.frame
    names(fetaldf)=c("chr", "start", "end", "State")
    fetaldf$Value = 1
    fetaldf %>% select(State, Value) -> fetaldf
    if (f == 1){
        fetalDF = fetaldf
    }else{
        fetalDF=rbind(fetalDF, fetaldf)
    }
}
fetalDF %>% group_by(State) %>% summarize(Val = sum(Value)) %>% filter(Val >= 2) -> fetalDF


adultBed=list.files("DiffStates/", pattern=paste0("^", tissue, "(.)*.", state, ".bed"))
for (f in 1:length(adultBed)){
    adultdf = fread(paste0("DiffStates/", adultBed[f]), header=F) %>% as.data.frame
    names(adultdf)=c("chr", "start", "end", "State")
    adultdf$Value = 1
    adultdf %>% select(State, Value) -> adultdf
    if (f == 1){
        adultDF = adultdf
    }else{
        adultDF=rbind(adultDF, adultdf)
    }
}
adultDF %>% group_by(State) %>% summarize(Val = sum(Value)) %>% filter(Val >= 2) -> adultDF

x=list(Fetal=as.vector(fetalDF$State), Adult=as.vector(adultDF$State))
venn.diagram(x,
            category.names = c(paste0("Fetal", tissue) , tissue),
            filename = paste0("DiffStates/", tissue, '.fetal_vs_adult.', state,'.tif'),
            imagetype = "tiff",
            height = 2400,
            width = 2400,
            resolution = 300,
            compression = "lzw",
            fill = c("#D6604D", '#4393C3'),
            sigdigs=1,
            fontfamily = "arial", cat.fontfamily="arial",
            cex = 2, cat.cex=2,
            margin=0.1
            )
fetalSpec = as.data.frame(setdiff(as.vector(fetalDF$State), as.vector(adultDF$State)))
adulSpec = as.data.frame(setdiff(as.vector(adultDF$State), as.vector(fetalDF$State)))
fwrite(fetalSpec, paste0("DiffStates/", tissue,".", state,".Fetal_specific.bed"), col.names=F, sep="\t", quote=F, row.names=F)
fwrite(adulSpec, paste0("DiffStates/", tissue,".", state,".Adult_specific.bed"), col.names=F, sep="\t", quote=F, row.names=F)
