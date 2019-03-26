library(compiler)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomicRanges)
#' findNDRs, based on Aaron Statham original code, plus Repitools bits (no R >3.5 version available)
#' Find NDRs in dSMF-seq data using a sliding window chi-squared test versus whole genome background
#'
#' @param x \code{GRanges} of methylation count data at GCH sites
#' @param samples \code{data.frame} describing the samples to discover NDRs for
#' @param p.cutoff Minimum -log10 p.value for a window to be considered significant
#' @param windowWidth Size of each windoe
#' @param windowBy Distance between the starts of consecutive windows
#' @param minSize Minimum size of an NDR when overlapping significant windows are pooled
#' @return \code{GRangesList} of NDRs found in each sample with the mean pvalue for each region
#' 
#' @importFrom compiler cmpfun
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges findOverlaps GRangesList values reduce width
#'

setGeneric("genomeBlocks", function(genome, ...) {standardGeneric("genomeBlocks")})
setMethod("genomeBlocks", "numeric",
          function(genome, chrs = names(genome), width = NULL, spacing = width)
          {
            if(is.null(width))
              stop("Block width not given.")
            
            chr.windows <- lapply(chrs, function(x) {
              centres <- seq.int(min(spacing/2,genome[x]),genome[x],spacing)
              GenomicRanges::GRanges(seqnames = names(genome[x]),
                      ranges = IRanges::IRanges(start = pmax(centres-width/2+1,1), 
                                       end=pmin(centres+width/2, genome[x])))
            })
            
            suppressWarnings(do.call(c, chr.windows))
          })

overlapSums <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(class(z)=="numeric" | class(z)=="integer")
  stopifnot(length(x)==length(z))
  ov <- as.matrix(GenomicRanges::findOverlaps(y, x))
  ovSums <- rep(as.numeric(NA), length(y))
  ovSums[unique(ov[,1])] <- IRanges::viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovSums
}
overlapMeans <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(class(z)=="numeric" | class(z)=="integer")
  stopifnot(length(x)==length(z))
  ov <- as.matrix(GenomicRanges::findOverlaps(y, x))
  ovMeans <- rep(as.numeric(NA), length(y))
  ovMeans[unique(ov[,1])] <- IRanges::viewMeans(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovMeans
}
fast.chisq <- compiler::cmpfun(function(x, p) {n <- rowSums(x)
  E <- cbind(n * p[1], n * p[2])
  STATISTIC <- rowSums((x - E)^2/E)
  PARAMETER <- 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  return(PVAL)
})


# Input data (x) is a genomic range for all Cs, with 
# total number of counts (or coverage, T column of the metadata)
# methylation counts(M column of the metadata)
# Samples names are taken from the data
# Order matters (for each sample, should have _T(total) then _M(methylated))

# Outputs a bigwig file with chi-square p-values versus average whole genome methylation
# Calls the "Nucleosome Depleted Regions" based on the set p.cutoff 
# Outputs a bifwig file with the NDRs

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

#Extract samples names from the data
column_names <- names(mcols(x))
M_columns <- grep("*_M", column_names, value=TRUE)
M_columns <- M_columns[order(M_columns)]
T_columns <- grep("*_T", column_names, value=TRUE)
T_columns <- T_columns[order(T_columns)]
stopifnot(length(M_columns)==length(T_columns))
samples_names <- substr(M_columns,1,nchar(M_columns)-2)
stopifnot(substr(T_columns,1,nchar(T_columns)-2)==samples_names)
samples <- data.frame(samples_names, T_columns, M_columns)
names(samples) <- c("Sample","T","M")
samples


#Parameters
windowWidth= 100 # window for detection (bin)
windowBy=20      # shift between successive windows
minSize=140      # minimal size of the detected NDR
p.cutoff=15      # cutoff for the chi-square p.value (10^-p.cutoff)

#x must have seqlengths defined
stopifnot(all(!is.na(seqlengths(x))))
# Define windows for pvalue calculation
windows <- genomeBlocks(GenomeInfoDb::seqlengths(x), width=windowWidth, spacing=windowBy)
GenomeInfoDb::seqlengths(windows) <- GenomeInfoDb::seqlengths(x)

for (j in 1:nrow(samples))
 {message(paste0("Processing ", samples$Sample[j]))
  # Counts of Ms and Cs in each window
  message(" - Counting mCs")
  windows.counts <- cbind("M"=overlapSums(x, windows, values(x)[[samples$M[j]]]))
  message(" - Counting Cs")
  windows.counts <- cbind(windows.counts, "T"=overlapSums(x, windows, values(x)[[samples$T[j]]])-
                            windows.counts[,"M"])
  
  # Genomic background ratio to compare to
  Csum <- sum(values(x)[[samples$M[j]]])
  Tsum <- sum(values(x)[[samples$T[j]]])-sum(values(x)[[samples$M[j]]])
  C.T <- c(Csum,Tsum)/(Tsum+Csum)

  # calculate chisq pvalue
  message(" - Calculating chisq pvalues")
  windows.pvals <- -log10(fast.chisq(windows.counts, C.T))
  windows.pvals[windows.pvals==Inf] <- max(windows.pvals[!windows.pvals==Inf])
  windows.pvals[is.na(windows.pvals)==TRUE] <- 0
  windows$score <- windows.pvals
  suppressWarnings(do.call(c, chr.windows))
  windows.for.bw <- GenomicRanges::disjoin(GenomicRanges::trim(GenomicRanges::resize(GenomicRanges::shift(windows,40),20)),with.revmap=TRUE)
  revmap <- mcols(windows.for.bw)$revmap
  r_scores <- extractList(mcols(windows)$score, revmap)
  mcols(windows.for.bw)<-NULL
  mcols(windows.for.bw)$score <- mean(r_scores)
  rtracklayer::export.bw(windows.for.bw, paste(samples$Sample[j],"_pvalues.bw", sep=""))
  # find regions that meet the cutoff on pvalue and size
  message(" - Finding significant regions")
  windows.cutoff <- reduce(windows[which(windows.pvals>p.cutoff)])
  windows.cutoff$score <- overlapMeans(windows[!is.na(windows.pvals)], windows.cutoff, windows.pvals[!is.na(windows.pvals)])
  #        windows.cutoff$p.mean <- overlapMeans(windows, windows.cutoff, windows.pvals, na.rm=TRUE)
  windows.cutoff <- windows.cutoff[width(windows.cutoff)>=minSize]
  rtracklayer::export.bw(windows.cutoff,paste(samples$Sample[j],"_NDRs_",p.cutoff,".bw",sep=""))
}

