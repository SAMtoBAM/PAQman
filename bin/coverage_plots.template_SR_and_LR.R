

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
###### SHORT-READ ALIGNMENT ###### 
###### ###### ###### ###### ###### 

#### read in the normalised coverage plots
SRcoverage=read.csv(file="PATHTOSRCOVERAGE", header=T, sep='\t' )
#### generate the plot of the relative coverage
## added lines donoting a 1X relative coverage and 0.5; no higher as I don't know if it'll be within the scale)
## can easily add higher by copying the geom_hline and changing the yintercept value
SRplot=ggplot(data=SRcoverage, aes(x=((start+end)/2), y=coverage_norm ))+
  geom_hline(yintercept = 1, linetype="dashed", colour="red")+
  geom_line()+
  facet_wrap(.~ contig , ncol=2, scales="free_y", )+
  theme_pubr()+
  scale_y_continuous(breaks = function(z) seq(0, range(z)[2], by = 0.5))+
  xlab("contig coordinates (bp)")+
  ylab("Genome-wide median normalised coverage")

#### output the plot as an svg
##can have issues with the plot size depending on the number of contig so calculate a plot height relative to that number
mmhigh=(length(unique(SRcoverage$contig))/2*50)
##output with width of just a bit less than A4 and DPI of 320
ggsave(filename = "PATHTOOUTPUT.SR.svg", plot = SRplot,
       width = 200, height = mmhigh, units = "mm", dpi = "retina")



###### ###### ###### ###### ###### 
###### LONG-READ ALIGNMENT  ###### 
###### ###### ###### ###### ###### 

#### read in the normalised coverage plots
LRcoverage=read.csv(file="PATHTOLRCOVERAGE", header=T, sep='\t' )
#### generate the plot of the relative coverage
## added lines donoting a 1X relative coverage and 0.5; no higher as I don't know if it'll be within the scale)
## can easily add higher by copying the geom_hline and changing the yintercept value
LRplot=ggplot(data=LRcoverage, aes(x=((start+end)/2), y=coverage_norm ))+
  geom_hline(yintercept = 1, linetype="dashed", colour="red")+
  geom_line()+
  facet_wrap(.~ contig , ncol=2, scales="free_y", )+
  theme_pubr()+
  scale_y_continuous(breaks = function(z) seq(0, range(z)[2], by = 0.5))+
  xlab("contig coordinates (bp)")+
  ylab("Genome-wide median normalised coverage")

#### output the plot as an svg
##can have issues with the plot size depending on the number of contig so calculate a plot height relative to that number
mmhigh=(length(unique(LRcoverage$contig))/2*50)
##output with width of just a bit less than A4 and DPI of 320
ggsave(filename = "PATHTOOUTPUT.LR.svg", plot = LRplot,
       width = 200, height = mmhigh, units = "mm", dpi = "retina")



###### ###### ###### ###### ###### ######
#### SHORT- AND LONG-READ ALIGNMENT #####
###### ###### ###### ###### ###### ######

## plot and save the relative coverage of the short and long reads

SRLRplot=ggplot()+
  geom_hline(yintercept = 1, linetype="dashed", colour="red")+
  geom_line(data=LRcoverage, aes(x=((start+end)/2), y=coverage_norm , colour="Long-read"))+
  geom_line(data=SRcoverage, aes(x=((start+end)/2), y=coverage_norm , colour="Short-read"))+
  facet_wrap(.~ contig , ncol=2, scales="free_y", )+
  theme_pubr()+
  scale_y_continuous(breaks = function(z) seq(0, range(z)[2], by = 0.5))+
  scale_color_simpsons()+
  xlab("contig coordinates (bp)")+
  ylab("Genome-wide median normalised coverage")

ggsave(filename = "PATHTOOUTPUT.SR_and_LR.svg", plot = SRLRplot,
       width = 200, height = mmhigh, units = "mm", dpi = "retina")