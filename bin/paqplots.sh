#!/bin/bash
set -euo pipefail

version="v1"

##paqman environment
#mamba create -n paqman bioconda::busco bioconda::merqury bioconda::quast bioconda::filtlong bioconda::seqtk bioconda::craq bioconda::seqkit conda-forge::r-gggenomes conda-forge::r-ggpubr conda-forge::r-ggsci conda-forge::r-svglite conda-forge::r-fmsb conda-forge::r-ggsci conda-forge::r-reshape2

##maybe to add mummer4


##to create the conda env using these tools and paqman
#conda env paqman > paqman.yml


paqplotspath=$( which comparison_paqman_plots.template.R )


##############################################################
##################### SETTING VARIABLES ######################
##############################################################

#default values, unless denoted when running PAQman
summary=""
list=""

prefix="paqplot"
output="paqplot_output"
help="nohelp"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-s|--summary)
	summary="$2"
	shift
	shift
	;;
	-l|--list)
	list="$2"
	shift
	shift
	;;
	-p|--prefix)
	prefix="$2"
	shift
	shift
	;;
	-o|--output)
	output="$2"
	shift
	shift
	;;
	-h|--help)
	echo "
	PAQman (version: ${version})
 
	paqplot -s summary_file.tsv
	OR
	paqplot -l list_of_summary_files.txt

	
	Required inputs:
	-s | --summary		A PAQman summary file with multiple assemblies combined and the same header
	-l | --list		A list of paths to multiple summary files
	
	Optional parameters:
	-p | --prefix		Prefix for output (default: paqplot)
	-o | --output		Name of output folder for all results (default: paqplot_output)
	-h | --help		Print this help message
	"
	exit
	;;
	esac
done

#creates error message and exits if these values are not/incorrectly assigned 
[[ $summary == "" && $list == "" ]] && echo "ERROR: Neither a path to a summary file or list of summary files was given. Provide one with -s or -l respectively." && exit

#some output to notify the use of some parameters
[[ $summary != "" ]] && echo "PAQman summary file given, using this file"
##set a full path for the summary file data to be read into the R script
[[ $summary != "" ]] && summarypath=$( realpath ${summary} )


##############################################################
#################### BEGINNING PAQman-plots ####################
##############################################################
echo "###################################################################"
echo "################## PAQman: Starting PAQman-plots ##################"
echo "###################################################################"
echo "################## PAQman: Step 1a: Collecting stats from summary files"


## assembly evaluations
mkdir ${output}


##create combined summary file if provided (-l) and -s was empty

if [[ ${summary} == "" && ${list} != "" ]]
then
echo "prefix;assembly;quast_#contigs;quast_#contigs>10kb;quast_assembly_size;quast_assembly_N50;quast_assembly_N90;quast_largest_contig;BUSCO_db;BUSCO_total;BUSCO_complete;BUSCO_complete_single;BUSCO_fragmented;BUSCO_missing;merqury_kmer_completeness(%);merqury_qv(phred);CRAQ_R-AQI(%);CRAQ_S-AQI(%);coverage_normal(%);telomeric_ends;telomeric_ends(%);T2T_contigs" | tr ';' '\t' > ${output}/combined.summary_stats.tsv
cat ${list} | while read file
do
tail -n+2 ${file} >> ${output}/combined.summary_stats.tsv
done
summarypath=$( realpath ${output}/combined.summary_stats.tsv )
fi

cd ${output}

##check if the combination of sample names and assembly names (if they are not just the same thing) creates unique names per sample
##first get an indicator if all the sample and assembly names are the same (0=same, 1=notsame)
allsame=$( tail -n+2 ${summarypath} | awk -F'\t' '$1!=$2 {exit 1}' ; echo $? )
##if the sample and assembly names are the same check just for uniqueness using this
if [[ "$allsame" -eq 0 ]]
then
dups=$( awk -F'\t' '{print $1}' "${summarypath}" | sort | uniq -d )
if [[ -n "$dups" ]]
then
echo "ERROR: Duplicate sample identifiers found: ${dups}"
echo "ERROR: Please use unique sample or assembly names (first two columns in summary files) to allow comparisons" && exit
fi
else
##if the sample and assembly names are NOT the same check just for uniqueness using a combination split by '-'
dups=$( awk -F'\t' '{print $1"-"$2}' "${summarypath}" | sort | uniq -d )
if [[ -n "$dups" ]]; then
echo "ERROR: Duplicate combined (sample-assembly) identifiers found: ${dups}" 
echo "ERROR: Please use unique sample or assembly names (first two columns in summary files) to allow comparisons" && exit
fi
fi


echo "################## PAQman: Step 1b: Plotting comparisons"
##covert paths in R script to those relevant for the coverage outputs generated here and then export the plots 
cat ${paqplotspath} | sed "s|PATHTOSUMMARY|${summarypath}|" | sed "s|PATHTOOUTPUT|./${prefix}|" > ./${prefix}.paqman_plots.R
Rscript ./${prefix}.paqman_plots.R


echo "################## PAQman: PAQman plots can be found here ${output}/*.svg"
echo "################## PAQman: All complete; thanks for using PAQman"

