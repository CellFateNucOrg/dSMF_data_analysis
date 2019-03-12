library(compiler)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(GenomicRanges)
#' findNDRs, based on Aaron Statham original code, plus Repitools bits (no R >3.5 version available)
#'#' Find NDRs in NOMe-seq data using a sliding window chi-squared test versus whole genome background
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
#' @importFrom Repitools genomeBlocks
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges findOverlaps GRangesList values reduce width
#'
#' @export
#'
#' @author Aaron Statham <a.statham@@garvan.org.au>
#' 

setGeneric("genomeBlocks", function(genome, ...) {standardGeneric("genomeBlocks")})
setMethod("genomeBlocks", "numeric",
          function(genome, chrs = names(genome), width = NULL, spacing = width)
          {
            if(is.null(width))
              stop("Block width not given.")
            
            chr.windows <- lapply(chrs, function(x) {
              centres <- seq.int(min(spacing/2,genome[x]),genome[x],spacing)
              GRanges(seqnames = names(genome[x]),
                      ranges = IRanges(start = pmax(centres-width/2+1,1), 
                                       end=pmin(centres+width/2, genome[x])))
            })
            
            suppressWarnings(do.call(c, chr.windows))
          })

overlapSums <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(class(z)=="numeric" | class(z)=="integer")
  stopifnot(length(x)==length(z))
  ov <- as.matrix(findOverlaps(y, x))
  ovSums <- rep(as.numeric(NA), length(y))
  ovSums[unique(ov[,1])] <- viewSums(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovSums
}
overlapMeans <- function(x, y, z, na.rm=FALSE) {
  stopifnot(class(x)=="GRanges")
  stopifnot(class(y)=="GRanges")
  stopifnot(class(z)=="numeric" | class(z)=="integer")
  stopifnot(length(x)==length(z))
  ov <- as.matrix(findOverlaps(y, x))
  ovMeans <- rep(as.numeric(NA), length(y))
  ovMeans[unique(ov[,1])] <- viewMeans(Views(z[ov[,2]], ranges(Rle(ov[,1]))), na.rm=na.rm)
  ovMeans
}
fast.chisq <- compiler::cmpfun(function(x, p) {n <- rowSums(x)
  E <- cbind(n * p[1], n * p[2])
  STATISTIC <- rowSums((x - E)^2/E)
  PARAMETER <- 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  return(PVAL)
})


setwd("C:/Users/pmeister/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/pmeister/Bolaji_dSMF/NDR_calls/")
#x <- readRDS("Amplicon_raw_methylation.rds")

x <- readRDS("dSMFproj_allCs_gw_counts.rds")
#Methylation is called and raw reads are now in 2 columns: 
#T (total, methylated and non methylated (=.cov))
#M (methylated, non converted (.C))
# Input data (x) is a genomic range for all C in a CG or GC context, with 
# total number of counts (or coverage, T column of the metadata)
# methylation counts(M column of the metadata)
# Order matters (for each sample, should have _T(total) then _M(methylated))


#Extract samples names from the data
column_names <- names(mcols(x))
column_names
samples_names <- substr(names(mcols(x)),1,nchar(names(mcols(x)))-2)[seq(1,length(names(mcols(x))),2)]
samples_names

#Make addition for mean 
mcols(x)[,5]<- mcols(x)[,1]+mcols(x)[,3]
mcols(x)[,6]<- mcols(x)[,2]+mcols(x)[,4]
samples_names <- c(samples_names,"N2_dSMFgw_16_20avg") 
names(mcols(x)) <- c(column_names, "N2_dSMFgw_16_20avg_T", "N2_dSMFgw_16_20avg_M")

samples <- data.frame(samples_names, names(mcols(x))[seq(1,length(names(mcols(x))),2)], names(mcols(x))[seq(2,length(names(mcols(x))),2)])
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
windows <- genomeBlocks(seqlengths(x), width=windowWidth, spacing=windowBy)

NDRs <- GRangesList(lapply(1:nrow(samples)
                           , function(j) {
  message(paste0("Processing ", samples$Sample[j]))
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
##TO DO Create a wiggle of the p-value on the windows for display
  # find regions that meet the cutoff on pvalue and size
  message(" - Finding significant regions")
  windows.cutoff <- reduce(windows[which(windows.pvals>p.cutoff)])
  windows.cutoff$p.mean <- overlapMeans(windows[!is.na(windows.pvals)], windows.cutoff, windows.pvals[!is.na(windows.pvals)])
  #        windows.cutoff$p.mean <- overlapMeans(windows, windows.cutoff, windows.pvals, na.rm=TRUE)
  windows.cutoff <- windows.cutoff[width(windows.cutoff)>=minSize]
}
))

NDRs_avg <- NDRs[[3]]
names(mcols(NDRs_avg))<- "score"
seqlengths(NDRs_avg) <- seqlengths(x)
export.bw(NDRs_avg,"N2_dSMF_NDR_16_20_avg.bw")
