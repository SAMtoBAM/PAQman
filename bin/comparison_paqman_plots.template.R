
################### ################### ################### ################### 
########################### READ IN LIBRARIES ################################
################### ################### ################### ################### 


#packages = c("ggplot2", "reshape2", "ggsci", "ggpubr", "devtools", "stringr")

##check if a package is installed, if so load, if not install
#package.check <- lapply(
#  packages,
#  FUN = function(x) {
#    if (!require(x, character.only = TRUE)) {
#      install.packages(x, dependencies = TRUE)
#      library(x, character.only = TRUE)
#    }
#  }
#)


##again for ggradar as it need devtools
#packages2=c("ggradar")
#package.check <- lapply(
#  packages2,
#  FUN = function(x) {
#    if (!require(x, character.only = TRUE)) {
#      devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
#      library(x, character.only = TRUE)
#    }
#  }
#)

library(ggplot2)
library(reshape2)
library(ggsci)
library(ggpubr)
library(stringr)
library(ggradar)


################### ################### ################### ################### 
###################### COMPARING SUMMARY STATISTICS ###########################
################### ################### ################### ################### 

### read in summary stats data
comparisons=read.csv(file="PATHTOSUMMARY", sep='\t', header=T, check.names=FALSE)

##place a combination of the strain and assembly name as the row name
#rownames(comparisons)=paste(comparisons$strain, comparisons$assembly, sep = "-")

##calculcate BUSCO complete counts as percentages
comparisons$"BUSCO_complete(%)"=(comparisons$BUSCO_complete/comparisons$BUSCO_total)*100
comparisons$"BUSCO_complete_single(%)"=(comparisons$BUSCO_complete_single/comparisons$BUSCO_total)*100


##set new variable as a combination of the strain and assembly name
comparisons$label=paste(comparisons$prefix, comparisons$assembly, sep = "-")

##order the label factor (alphanumerically) for each plot
list=str_sort(comparisons$label)
comparisons$label=factor(comparisons$label, levels=list)

##plots of raw values
## will generate a set of plots, one being a radar plot of percentages, the rest being the absolute values in side by side columns

##select only the variables that will be compared
##here we have
## quast_#contigs
## quast_assembly_size
## quast_assembly_N50
## BUSCO_complete
## merqury_completeness(%)
## merqury_qv(phred)
## CRAQ_average_CRE(%)
## CRAQ_average_CSE(%)
## coverage_normal(%)
## telomeric_ends
## t2t_contigs
comparisonsrad=subset(comparisons, select = c("BUSCO_complete(%)", "BUSCO_complete_single(%)", "merqury_completeness(%)", "CRAQ_average_CRE(%)", "CRAQ_average_CSE(%)", "coverage normal(%)", "telomeric_ends(%)"))
#comparisonsrad=comparisons[ , c(25, 24, 23, 15, 17, 18, 19, 21), FALSE]
##tidy up the headers
names(comparisonsrad) = gsub(pattern = "merqury_", replacement = "", x = names(comparisonsrad))
names(comparisonsrad) = gsub(pattern = "quast_", replacement = "", x = names(comparisonsrad))
names(comparisonsrad) = gsub(pattern = "CRAQ_average_", replacement = "", x = names(comparisonsrad))

absplot1=ggradar(comparisonsrad, axis.label.size = 3, legend.text.size = 6, legend.position = "left", group.point.size = 3 , group.line.width = 1, gridline.mid.colour = "grey", grid.label.size = 5, background.circle.colour = "grey90", gridline.mid.linetype = 8 , gridline.max.linetype = 8, grid.max = 100.1, grid.mid=50,  values.radar = c("", "50%", ""))+scale_color_aaas()+coord_cartesian(clip = "off")+theme(plot.margin = margin(0, 5, 0, 5, 'cm'))

##now set up for the lollipop plots
comparisonsdot=subset(comparisons, select = c("quast_#contigs", "quast_#contigs>10kb", "quast_assembly_size", "quast_assembly_N50", "quast_assembly_N90", "quast_largest_contig", "merqury_qv(phred)", "telomeric_ends", "t2t_contigs"))
#comparisonsdot=comparisons[ , c(25, 3, 4, 5, 6, 7, 8, 16, 20, 22), FALSE]
##remove some of the naming conventions as the labels are too big with them
names(comparisonsdot) = gsub(pattern = "merqury_", replacement = "", x = names(comparisonsdot))
names(comparisonsdot) = gsub(pattern = "quast_", replacement = "", x = names(comparisonsdot))

##melt the dataframe
comparisonsdot2=melt(comparisonsdot)

##plot 
absplot2=ggdotchart(comparisonsdot2, x = "label", y="value", sorting = "none", rotate=T, color="label", facet.by="variable", scales="free_x", nrow=1, add="segment", legend = "none") + theme_pubr(x.text.angle = 45, legend = "none") + xlab("Assembly")+scale_color_aaas()+theme(axis.text.y = element_text(size=10))


##combine the two plots together
absplot3=ggarrange(absplot1, absplot2, ncol=1, heights = c(3,2))

##save the plot
ggsave(filename = "PATHTOOUTPUT.raw_values.svg", plot = absplot3,
       width = 500, height = 250, units = "mm", dpi = "retina")



####RELATIVE VALUES PAQPLOT

##select only the variables that will be compared (all values will be made relative regardless)
##here we have
## quast_#contigs
## quast_assembly_size
## quast_assembly_N50
## BUSCO_complete(%)
## BUSCO_complete_single(%)
## merqury_completeness(%)
## merqury_qv(phred)
## CRAQ_average_CRE(%)
## CRAQ_average_CSE(%)
## coverage_normal(%)
## telomeric_ends
## t2t_contigs
comparisonstemp1=subset(comparisons, select = c("quast_assembly_size", "quast_assembly_N50", "quast_#contigs", "merqury_qv(phred)", "merqury_completeness(%)", "CRAQ_average_CRE(%)", "CRAQ_average_CSE(%)", "coverage normal(%)", "telomeric_ends", "t2t_contigs", "BUSCO_complete(%)", "BUSCO_complete_single(%)"))
#comparisonstemp1=comparisons[ , c(5, 6, 3, 15, 16, 17, 18, 19, 20, 22, 23, 24), FALSE]
##remove some of the naming conventions as the labels are too big with them
names(comparisonstemp1) = gsub(pattern = "merqury_", replacement = "", x = names(comparisonstemp1))
names(comparisonstemp1) = gsub(pattern = "quast_", replacement = "", x = names(comparisonstemp1))
names(comparisonstemp1) = gsub(pattern = "CRAQ_average_", replacement = "", x = names(comparisonstemp1))


##get a label to be placed in the first column (taken by ggradar as the label)
label=as.data.frame(paste(comparisons$prefix, comparisons$assembly, sep = "-"))
##rename column name
colnames(label)="label"


##calculate the relative values of all samples for each column, dividing by the max for the column
comparisonstemp2=as.data.frame(apply(comparisonstemp1,2,function(x){x/max(x)}, simplify = F), check.names = F)

##combine the labels (in the first column) with the relative values
comparisonsrel=cbind(label, comparisonstemp2)
list=str_sort(comparisonsrel$label)
comparisonsrel$label=factor(comparisonsrel$label, levels=list)

##plot with ggradar
relplot=ggradar(comparisonsrel, axis.label.size = 3, legend.text.size = 6, legend.position = "top", group.point.size = 3 , group.line.width = 1, gridline.mid.colour = "grey", grid.label.size = 5, background.circle.colour = "grey90", gridline.mid.linetype = 8 , gridline.max.linetype = 8, values.radar = c("", "0.5", ""))+scale_color_aaas()+theme(legend.direction = "vertical")

ggsave(filename = "PATHTOOUTPUT.relative_values.svg", plot = relplot,
       width = 200, height = 200, units = "mm", dpi = "retina")
