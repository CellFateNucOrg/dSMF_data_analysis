library(Biostrings)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(stats)

path = "C:/Users/pmeister/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/pmeister/Bolaji_dSMF/NDR_calls/" 

setwd(path)
#x <- readRDS("Amplicon_raw_methylation.rds")
all_Cs <- readRDS("dSMFproj_allCs_gw_counts.rds")
genome=Celegans

#Remove Cs in a non CG or GC context
CGs<-Biostrings::vmatchPattern("CG", genome)
CGs<-GenomicRanges::reduce(CGs, ignore.strand=TRUE)
GCs<-Biostrings::vmatchPattern("GC", genome)
GCs<-GenomicRanges::reduce(GCs,ignore.strand=TRUE)
meCs<-GenomicRanges::union(CGs,GCs)
meCs<-GenomicRanges::reduce(meCs)
x <- all_Cs[subjectHits(GenomicRanges::findOverlaps(meCs, all_Cs,ignore.strand=TRUE))]

#Build Bsseq object

M <- as.matrix(alldat[,ix.X, drop=FALSE])
Cov <- as.matrix(alldat[,ix.N, drop=FALSE])
colnames(M) <- colnames(Cov) <- sampleNames


N2_dSMFv020gw <- BSseq(chr=as.vector(seqnames(x)), 
              pos=start(x),
              M=matrix(mcols(x)$N2_dSMFv020gw_M), 
              Cov=matrix(mcols(x)$N2_dSMFv020gw_T))
head(N2_dSMFv020gw, n=10)
sm50k_N2_dSMFv020gw<- BSmooth(N2_dSMFv020gw, h= 50000)
wig_for_export <- getBSseq(sm_N2_dSMFv020gw,type="gr")
seqlengths(wig_for_export)<- seqlengths(x)
wig_for_export$score <- getMeth(sm_N2_dSMFv020gw)
export.bw(wig_for_export,"sm_N2_dSMFv020gw_50k.bw")

mean(getBSseq(sm_N2_dSMFv020gw,type="coef"))


sm_N2_dSMFv020gw<- BSmooth(N2_dSMFv020gw, h= 40)
wig_for_export <- getBSseq(sm50_N2_dSMFv020gw,type="gr")
seqlengths(wig_for_export)<- seqlengths(x)
wig_for_export$score <- 1-getMeth(sm50_N2_dSMFv020gw)
export.bw(wig_for_export,"sm_N2_dSMFv020gw_100.bw")


sm40_N2_dSMFv020gw<- BSmooth(N2_dSMFv020gw, h= 40)
wig_for_export <- getBSseq(sm40_N2_dSMFv020gw,type="gr")
seqlengths(wig_for_export)<- seqlengths(x)
wig_for_export$score <- 1-getMeth(sm40_N2_dSMFv020gw)
export.bw(wig_for_export,"sm_N2_dSMFv020gw_40.bw")
