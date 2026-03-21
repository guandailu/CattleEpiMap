library(gkmSVM) 
library(BSgenome.Btaurus.UCSC.bosTau9.masked)
args <- commandArgs(trailingOnly = TRUE)
mark=args[1]
tissue=args[2]

### prepare positive and negtive sets

peaks=paste0("data/", mark, "/", tissue,"/idr_output/", mark,"_", tissue,"_overlap_peak.bed")
#peaks=paste0("../chrombpnet/ATAC_bams/", tissue, "/idr_output/", mark,"_", tissue,"_overlap_peak.bed")
syscmd=paste0("cat ", peaks, " | awk -v OFS=\"\t\" '{print\"chr\"$1, $2,$3}' > data/", mark, "/", tissue, "/", mark,"_", tissue, "_overlap_peak.bed")
#syscmd=paste0("cat ", peaks, " | awk -v OFS=\"\t\" '{print $1, $2,$3}' > data/", mark, "/", tissue, "/", mark,"_", tissue, "_overlap_peak.bed")
system(syscmd)


genNullSeqs(inputBedFN=paste0("data/", mark,"/", tissue,"/", mark,"_", tissue,"_overlap_peak.bed"), 
            nMaxTrials=10, 
            xfold=1, 
            genome=BSgenome.Btaurus.UCSC.bosTau9.masked, 
            outputPosFastaFN=paste0("output/", mark,"/", tissue,"/", mark,"_", tissue,"_overlap_peak.pos1x.fa"), 
            outputBedFN=paste0("output/", mark,"/", tissue,"/", mark,"_", tissue,"_overlap_peak.neg1x.bed"), 
            outputNegFastaFN=paste0("output/", mark,"/", tissue,"/", mark,"_", tissue,"_overlap_peak.neg1x.fa"))

