#path = "C:/Users/pmeister/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/pmeister/dSMF_data_analysis/" 

#setwd(path)
#x <- readRDS("Amplicon_raw_methylation.rds")
#data <- readRDS("sm_DSMF_hotsitesA.RDS")
library(ggplot2)
library(genomation)
library(tidyr)
####################################################
# Plotting function for genomation ScoreMatrix     #
####################################################
#' @importFrom ggplot2 ggplot
#' @importFrom genomation ScoreMatrix
#' @importFrom tidyr gather
#' @param data   \code {genomation} ScoreMatrix of dSMF values[0-1]
#' @param myYlab {text} Y label for graph
#' @param myXlab {text} X label for graph
#' @param feature_label {text} middle value on which windowns are aligned to
#' @param title {text} graph titleDistance between the starts of consecutive windows
#' 
plotAveragedSMF<-function(data,         #genomation ScoreMatrix
                          myYlab,       #label for Y axis (type of loci)
                          myXlab="CpG/GpC position", #label for X axis
                          feature_label="TSS",
                          title=NULL)
  {
#Transform genomation matrix to data frame
df <- as.data.frame(data)
#Remove 0 values
df[df==0] <- NA
#rename columns, get number of columns
colnames(df)<-as.numeric(gsub("V","",colnames(df)))
width_data <- max(as.numeric(colnames(df)))
#name samples
rownames(df)<- as.numeric(rownames(df))
#prepare data for plotting
d<- tidyr::gather(df,key=position, value=methylation)
d$molecules<-rownames(df)
d$position<-as.numeric(d$position)

p <- ggplot2::ggplot(d,aes(x=position,y=molecules,width=1)) +
  ggplot2::geom_tile(aes(width=1,fill=methylation)) +
  ggplot2::scale_fill_gradientn(values = c(0,1),
                       colors=c("yellow", "red"),
                       na.value="black") +
  ggplot2::theme(panel.background = element_rect(fill="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =  element_blank(),
        plot.title = element_text(face = "bold",hjust = 0.5)) +
 # ggplot2::guides(colour=guide_colorbar(),
  #       fill = guide_legend(title="dSMF", 
   #                          title.position="top")) + 
  ggplot2::scale_x_continuous(name=myXlab, 
                     limits = c(0,width_data),
                     expand=c(0,0),
                     breaks=c(0,width_data/2,width_data),
                     labels=c(paste0("-",width_data/2),feature_label,paste0(width_data))) + 
  ggplot2::ylab(myYlab) + 
  ggplot2::ggtitle(title)
return(p)
}
