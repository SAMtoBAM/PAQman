

###### LIBRARIES ###### 
##for plotting
library(ggplot2)
library(ggpubr)
library(ggsci)
##to save plot as svg
library(svglite)


## This script is to automate the plotting of read coverage using short- and long-reads seperately
## Just need to change the strings in place of the path names
## The values for coverage are the median value of a 30kb (10kb sliding) bin to reduce the file size and smooth the alignment
## For plotting only the normalised coverage is used. This was calculcated by dividing the actual coverage by the genome wide median


###### ###### ###### ###### ###### 
###### LONG-READ ALIGNMENT  ###### 
###### ###### ###### ###### ###### 

#### read in the normalised coverage plots
LRcoverage=read.csv(file="PATHTOLRCOVERAGE", header=T, sep='\t' )

##make sure that each contig has at least >1 data point
LRcoverage_filtered <- subset(LRcoverage, ave(seq_along(contig), contig, FUN = length) > 1)

#### generate the plot of the relative coverage
## added lines donoting a 1X relative coverage and 0.5; no higher as I don't know if it'll be within the scale)
## can easily add higher by copying the geom_hline and changing the yintercept value
LRplot=ggplot(data=LRcoverage_filtered, aes(x=((start+end)/2), y=coverage_norm ))+
  geom_hline(yintercept = 1, linetype="dashed", colour="red")+
  geom_line()+
  facet_wrap(.~ contig , ncol=2, scales="free_y", )+
  theme_pubr()+
  scale_y_continuous(breaks = function(z) seq(0, range(z)[2], by = 0.5))+
  xlab("contig coordinates (bp)")+
  ylab("Genome-wide median normalised coverage")

#### output the plot as an svg
##can have issues with the plot size depending on the number of contig so calculate a plot height relative to that number
mmhigh=(length(unique(LRcoverage_filtered$contig))/2*50)
##output with width of just a bit less than A4 and DPI of 320
ggsave(filename = "PATHTOOUTPUT.LR.svg", plot = LRplot,
       width = 200, height = mmhigh, units = "mm", dpi = "retina", limitsize = F)

