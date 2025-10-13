#!/bin/bash
set -euo pipefail

version="v1.0.1"

LRcoverageRpath=$( which coverage_plots.template_LR.R )
LRSRcoverageRpath=$( which coverage_plots.template_SR_and_LR.R )


##############################################################
##################### SETTING VARIABLES ######################
##############################################################

#default values, unless denoted when running PAQman
assembly=""
longreads=""
platform="ont"
pair1=""
pair2=""
shortreads="no"
buscodb="eukaryota"
threads="1"
window="30000"
slide="10000"
telomererepeat="TTAGGG"
prefix="paqman"
output="paqman_output"
sequences="contigs"
LRmeryldb=""
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-a|--assembly)
	assembly="$2"
	shift
	shift
	;;
	-l|--longreads)
	longreads="$2"
	shift
	shift
	;;
	-x|--platform)
	platform="$2"
	shift
	shift
	;;
	-1|--pair1)
	pair1="$2"
	shift
	shift
	;;
	-2|--pair2)
	pair2="$2"
	shift
	shift
	;;
	-b|--buscodb)
	buscodb="$2"
	shift
	shift
	;;
	-t|--threads)
	threads="$2"
	shift
	shift
	;;
	-w|--window)
	window="$2"
	shift
	shift
	;;
	-s|--slide)
	slide="$2"
	shift
	shift
	;;
	-r|--telomererepeat)
	telomererepeat="$2"
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
	-seq|--sequences)
	sequences="$2"
	shift
	shift
	;;	
	-LRmdb|--LRmeryldb)
	LRmeryldb="$2"
	shift
	shift
	;;
	-c|--cleanup)
	cleanup="$2"
	shift
	shift
	;;
	-h|--help)
	echo "
	PAQman (version: ${version})
 
	paqman.sh -a assembly.fa -l long-reads.fq.gz
	
	Required inputs:
	-a | --assembly     Genome assemly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz)
	-l | --longreads    Long reads used for assembly in fastq or fasta format  (*.fa / *.fasta / *.fna / *.fastq / *.fq) and can be gzipped (*.gz)

	Recommended inputs:
	-x | --platform     Long-read technology to determine mapping mapping parameters. Choose between 'ont' or 'pacbio-hifi' or 'pacbio-clr' (default: ont)
	-b | --buscodb      Name of BUSCO database to be used (default: eukaryota)
	-t | --threads      Number of threads for tools that accept this option (default: 1)
	-r | --repeat       Telomeric repeat pattern (default: TTAGGG)
 	-1 | --pair1        Paired end illumina reads in fastq format; first pair. Used by Merqury, CRAQ and coverage analysis (Recommended). Can be gzipped (*.gz)
	-2 | --pair2        Paired end illumina reads in fastq format; second pair. Used by Merqury, CRAQ and coverage analysis (Recommended). Can be gzipped (*.gz)
	
	Optional parameters:
	-w | --window       Number of basepairs for window averaging for coverage (default: 30000)
	-s | --slide        Number of basepairs for the window to slide for coverage (default: 10000)
	-p | --prefix       Prefix for output (default: name of assembly file (-a) before the fasta suffix)
	-o | --output       Name of output folder for all results (default: paqman_output)
	-seq | --sequences	Whether or not to use scaffolds or contigs; provide 'scaffolds' to not break the assembly at N's (default: contigs)
	-LRmdb | --LRmeryldb	A precomputed meryl database for your long-read dataset. Will generate the meryl db if not provided
	-c | --cleanup      Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: yes)
	-h | --help         Print this help message
	"
	exit
	;;
	*)  # catch invalid args
    echo "ERROR: Unknown option: '$1'"
    echo "Run 'paqman.sh -h' to see valid options"
    exit 1
    ;;
	esac
done



##reset the prefix if not reassigned from 'paqman' the to the prefix of the assembly
[[ $prefix == "paqman" ]] && prefix=$( echo $assembly | awk -F "/" '{print $NF}' | sed 's/\.fasta\.gz$//' | sed 's/\.fa\.gz$//' | sed 's/\.fasta$//' | sed 's/\.fa$//' | sed 's/\.fna$//' )
##set that we do have shortreads available
[[ $pair1 != "" && $pair2 != ""  ]] && shortreads="yes"


#creates error message and exits if these values are not/incorrectly assigned 
[[ $assembly == "" ]] && echo "ERROR: Path to assembly not found, assign using -a" && exit
[[ $longreads == "" ]] && echo "ERROR: Path to long-reads not found, assign using -l" && exit
[[ $platform != "ont" && $platform != "pacbio-hifi" && $platform != "pacbio-clr" ]] && echo "ERROR: Incorrect platform assigned. Please select from 'ont', 'pacbio-hifi' or 'pacbio-clr' using -x" && exit

#some output to notify the use of some parameters
echo "Using ${platform} parameters for long-read alignment"
[[ $pair1 == "" || $pair2 == ""  ]] && echo "Missing short-reads as input; using only long-reads for assessment"
echo "Using ${buscodb} database for BUSCO assessment"

##assembly name without the suffix (can be the prefix too if not set)
assembly2=$( echo $assembly | awk -F "/" '{print $NF}' | sed 's/\.fasta\.gz$//' | sed 's/\.fa\.gz$//' | sed 's/\.fasta$//' | sed 's/\.fa$//' | sed 's/\.fna$//' )

##get paths for the reads and assembly
assemblypath=$( realpath ${assembly} )
longreadpath=$( realpath ${longreads} )
[[ $pair1 != "" || $pair2 != ""  ]] && pair1path=$( realpath ${pair1} ) && pair2path=$( realpath ${pair2} )

#check if the files given actually exist
[ ! -f "${assemblypath}" ] && echo "ERROR: Cannot find path to assembly file provided by -a; check path is correct and file exists" && exit
[ ! -f "${longreadpath}" ] && echo "ERROR: Cannot find path to long-reads provided by -l; check path is correct and file exists" && exit

if [[ $shortreads == "yes" ]]
then
[ ! -f "${pair1path}" ] && echo "ERROR: Cannot find path to short-reads provided by -1; check path is correct and file exists" && exit
[ ! -f "${pair2path}" ] && echo "ERROR: Cannot find path to short-reads provided by -2; check path is correct and file exists" && exit
fi

##check if meryl database for the long-reads was provided and if so generate a absolute path to be providede to Merqury
##exit if cannot find provided file
if [[ $LRmeryldb != "" ]]
then
LRmeryldbpath=$( realpath ${LRmeryldb} )
[ ! -f "${LRmeryldbpath}" ] && echo "ERROR: Cannot find path to meryl database provided by --LRmeryldb; check path is correct and file exists" && exit
fi


##check that the options are proper
[[ $sequences != "contigs" && $sequences != "scaffolds" ]] && echo "ERROR: --sequences option neither 'contigs' (default) or 'scaffolds'" && exit
[[ $platform != "ont" && $platform != "pacbio-hifi" && $platform != "pacbio-clr" ]] && echo "ERROR: --platform option needs to be 'ont' (default) or 'pacbio-hifi' or 'pacbio-clr'" && exit
[[ $cleanup != "yes" && $cleanup != "no" ]] && echo "ERROR: --cleanup option neither 'yes' (default) or 'no'" && exit

[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

##############################################################
#################### BEGINNING EVALUATION ####################
##############################################################
echo "#################################################"
echo "################ Starting PAQman ################"
echo "#################################################"
echo "$(date +%H:%M) ########## Step 1: Organising Input"

## assembly evaluations
mkdir ${output}
cd ${output}

##organise the assemblies
##if they are compressed then create a new version uncompressed
##if the sequences option is switched to scaffold then leave the sequences as is

## make symbolic links to the input files in the output folder
##unzip the assembly if it was compressed otherwise just create a symobolic link
if [[ ${sequences} == "scaffolds" ]]
then
echo "NOTE: Keeping sequences as is, not breaking scaffolds"
if [[ ${assembly} =~ ".gz"$ || ${assembly} =~ ".gzip"$ ]]
then
zcat ${assemblypath} | seqkit sort  > ${assembly2}.fa
assembly="${assembly2}.fa"
else
cat ${assemblypath} | seqkit sort  > ${assembly2}.fa
assembly="${assembly2}.fa"
fi
fi


##if the default option for sequences is to use contigs 
##if this default option is maintained then the assembly will be broken into contigs at all sets of Ns
##save the new assembly and use this for everything

if [[ ${sequences} == "contigs" ]]
then
echo "NOTE: Breaking any scaffolds into contigs for all evaluations"
if [[ ${assembly} =~ ".gz"$ || ${assembly} =~ ".gzip"$ ]]
then
zcat ${assemblypath} | awk -F " " '{print $1}' | seqkit fx2tab | awk -F'\t' '{
        header = $1;
        sequence = $2;
        split(sequence, contigs, /N+/);
        for (i=1; i<=length(contigs); i++) {
            if (length(contigs[i]) > 0) { # Ensure non-empty contigs
                print ">" header "_contig" i;
                print contigs[i];
            }
        }
    }' | seqkit sort  > ${assembly2}.fa
assembly="${assembly2}.fa"
else
cat ${assemblypath} | awk -F " " '{print $1}' | seqkit fx2tab | awk -F'\t' '{
        header = $1;
        sequence = $2;
        split(sequence, contigs, /N+/);
        for (i=1; i<=length(contigs); i++) {
            if (length(contigs[i]) > 0) { # Ensure non-empty contigs
                print ">" header "_contig" i;
                print contigs[i];
            }
        }
    }' | seqkit sort  > ${assembly2}.fa
assembly="${assembly2}.fa"
fi
fi


## make symbolic links to the input files in the output folder
ln -sf ${longreadpath} ./
[[ $shortreads == "yes" ]] && ln -sf ${pair1path} ./
[[ $shortreads == "yes" ]] && ln -sf ${pair2path} ./
## create modified variable for that calls directly the local symlink
longreads2=$( echo ${longreads} | awk -F "/" '{print $NF}' )
pair12=$( echo ${pair1} | awk -F "/" '{print $NF}' )
pair22=$( echo ${pair2} | awk -F "/" '{print $NF}' )



########################## QUAST ##########################
echo "$(date +%H:%M) ########## Step 2: Running Quast"
quast -o ./quast ${assembly} > quast.log
##move quast log to quast output folder
mv quast.log quast/
## we are just interested in the summary tsv file 'quast/report.tsv'


########################## BUSCO (using the ${buscodb} dataset) ##########################
echo "$(date +%H:%M) ########## Step 3: Running BUSCO"
busco -i ${assembly} -o ./busco  -l ${buscodb} --mode genome -c ${threads} >  busco.log
##move log to busco output folder
mv busco.log busco/

if [[ $cleanup == "yes" ]]
then
##remove output from busco that takes a lot of space
rm -r ./busco/run_${buscodb}*/*_output
rm -r ./busco/run_${buscodb}*/busco_sequences
rm -r ./busco/logs
[ -e "./busco/tmp" ] && rm -r ./busco/tmp
##remove the downloaded busco database
rm -r busco_downloads
fi
## we are just interested in the summary txt file 'busco/short_summary.*.txt'


########################## Merqury ##########################
echo "$(date +%H:%M) ########## Step 4a: Generating k-mer distribution"
if [[ $shortreads == "yes" ]]
then
echo "NOTE: Creating meryl k-mer database using short-reads"
## calculcate the kmer profile of the raw reads before assembly
## this profile will be compared to the resulting assembly to calculate completeness, i.e. how many of the good quality kmers are captured in the assembly
## here we can use JUST the illumina dataset and always compare to this dataset
meryl t=${threads} memory=15 k=18 count output reads.meryl ${pair12} ${pair22}
else 

if [[ $LRmeryldb != "" ]]
then
echo "NOTE: Using precomputed meryl k-mer database for long-reads" 
else
echo "NOTE: Creating meryl k-mer database using long-reads"
##same but instead using the long-read data due to an absence of short reads
meryl t=${threads} memory=15 k=18 count output reads.meryl ${longreads2} 
fi

fi

echo "$(date +%H:%M) ########## Step 5b: Running Merqury"

mkdir ./merqury
if [[ $LRmeryldb != "" ]]
then
merqury.sh ${LRmeryldbpath} ${assembly} ${prefix}.merqury
else
merqury.sh reads.meryl ${assembly} ${prefix}.merqury
fi
## just want to keep these two output files with important stats on error rate and completeness (respectively)
mv ${prefix}.merqury.qv ./merqury/
mv ${prefix}.merqury.completeness.stats ./merqury/
if [ -f "${prefix}.merqury.spectra-asm.fl.png" ] ; then mv ${prefix}.merqury.*.png ./merqury/ ; fi
if [ -f "${prefix}.merqury.${prefix}.spectra-cn.fl.png" ] ; then mv ${prefix}.merqury.${prefix}.*.png ./merqury/ ; fi
mv logs merqury/ 
if [ -f "${prefix}.merqury.spectra-asm.hist" ] ; then mv ${prefix}.merqury.*.hist ./merqury/ ; fi
if [ -f "${prefix}.merqury.${prefix}.qv" ] ; then mv ${prefix}.merqury.${prefix}.qv ./merqury/ ; fi


if [[ $cleanup == "yes" ]]
then
## remove the rest of the files we don't care about to clean up the directory
if [ -f "reads.filt" ] ; then rm reads.filt ; fi
if [ -f "reads.hist" ] ; then rm reads.hist ; fi
if [ -f "reads.hist.ploidy" ] ; then rm reads.hist.ploidy ; fi
if [ -f "${prefix}_only.bed" ] ; then rm ${prefix}_only.bed ; fi
if [ -f "${prefix}_only.wig" ] ; then rm ${prefix}_only.wig ; fi
rm -r *.meryl
fi



########################## CRAQ ##########################
echo "$(date +%H:%M) ########## Step 5a: Downsampling for 50X long-reads for CRAQ assessment"
## redownsample the dataset for just 50X of the longest as to remove shorter reads from confusing CRAQ
## get genome size based on input genome
genomesize=$( cat ./quast/report.tsv  | awk -F "\t" '{if(NR == 15) print $2}' )
target=$( echo $genomesize | awk '{print $1*50}' )
##now run filtlong with the settings
filtlong -t ${target} --length_weight 5 ${longreads2} | gzip > longreads.filtlong50x.fq.gz
##run craq
echo "$(date +%H:%M) ########## Step 5b: Running CRAQ"
if [[ $shortreads == "yes" ]]
then
[[ $platform == "ont" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -ngs ${pair12},${pair22} -x map-ont --thread ${threads} > craq.log
[[ $platform == "pacbio-hifi" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -ngs ${pair12},${pair22} -x map-hifi --thread ${threads} > craq.log
[[ $platform == "pacbio-clr" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -ngs ${pair12},${pair22} -x map-pb --thread ${threads} > craq.log
mv craq.log craq/
if [[ $cleanup == "yes" ]]
then
rm -r ./craq/SRout
fi
else
[[ $platform == "ont" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -x map-ont --thread ${threads} > craq.log
[[ $platform == "pacbio-hifi" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -x map-hifi --thread ${threads} > craq.log
[[ $platform == "pacbio-clr" ]] && craq -D ./craq -g ${assembly} -sms longreads.filtlong50x.fq.gz -x map-pb --thread ${threads} > craq.log
mv craq.log craq/
fi
## remove large intermediate files
if [[ $cleanup == "yes" ]]
then
rm -r ./craq/LRout
##remove read subset used for alignment
rm longreads.filtlong50x.fq.gz
fi
## from the output we are just interested in the summary file 'craq/runAQI_out/out_final.Report' for the summary stats (second line from the top is the genome wide average)
## for the position of structural errors, 'craq/runAQI_out/strER_out/out_final.CSE.bed'




########################## READ COVERAGE ##########################
echo "$(date +%H:%M) ########## Step 6a: Running whole genome read alignment"
## using the short-read data we can look at genome wide coverage on the assembly
## first index the assembly using the aligner
mkdir ./coverage
bwa index ${assembly}

## the window size for median-average binning (in bp)
window2=$( echo $window | awk '{print $0/1000}' )
## distance for the window to move before recalculating the window (in bp)
slide2=$( echo $slide | awk '{print $0/1000}' )
## create a bed file use for binning
samtools faidx ${assembly}
cat ${assembly}.fai | cut -f1-2 > ${assembly}.bed
## now split the reference bed file up into the above prescribed bins
bedtools makewindows -w ${window} -s ${slide} -g ${assembly}.bed > ./coverage/${prefix}.${window2}kbwindow_${slide2}kbslide.bed
## get the coverage file (slightly modify by giving a range for the single basepair coverage value) then use bedtools map to overlap with the reference derived window file to create median-averaged bins

###RUNNING THE LONG-READ ALIGNMENT
[[ $platform == "ont" ]] && minimap2 --secondary=no -ax map-ont -t ${threads} ${assembly} ${longreads2} | samtools sort -@ 4 -o ./coverage/${prefix}.minimap.sorted.bam -
[[ $platform == "pacbio-hifi" ]] && minimap2 --secondary=no -ax map-hifi -t ${threads} ${assembly} ${longreads2} | samtools sort -@ 4 -o ./coverage/${prefix}.minimap.sorted.bam -
[[ $platform == "pacbio-clr" ]] && minimap2 --secondary=no -ax map-pb -t ${threads} ${assembly} ${longreads2} | samtools sort -@ 4 -o ./coverage/${prefix}.minimap.sorted.bam -

## get the coverage
bedtools genomecov -d -split -ibam ./coverage/${prefix}.minimap.sorted.bam | gzip > ./coverage/${prefix}.minimap.sorted.cov.tsv.gz

## calculate the median across the whole genome using the exact basepair
medianLR=$( zcat ./coverage/${prefix}.minimap.sorted.cov.tsv.gz | awk '{if($3 != "0") print $3}' | sort -n | awk '{ a[i++]=$1} END{x=int((i+1)/2); if(x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' )
##calculate the binned median coverage and normalise each bin value by the genome wide median coverage
echo "contig;start;end;coverage_abs;coverage_norm" | tr ';' '\t' > ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv
zcat ./coverage/${prefix}.minimap.sorted.cov.tsv.gz  | awk '{print $1"\t"$2"\t"$2"\t"$3}' | bedtools sort |\
bedtools map -b - -a ./coverage/${prefix}.${window2}kbwindow_${slide2}kbslide.bed -c 4 -o median | awk -v median="$medianLR" '{print $0"\t"$4/median}' >> ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv

if [[ $cleanup == "yes" ]]
then
##remove the alignment file due to size
rm ./coverage/${prefix}.minimap.sorted.bam
##remove the coverage files looking at everybase pair due to size
rm ./coverage/${prefix}.minimap.sorted.cov.tsv.gz
##remove bed file generated from assembly
rm ${assembly}.bed
fi

## generate genome wide plots using the final file './coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.coverage_normalised.tsv'
##using a template 
LRpath=$( realpath ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv )

##only run this if short reads are provided
if [[ $shortreads == "yes" ]]
then
## align the short reads (filtering for only primary alignments -F 0x100 : removes secondary)
#bwa mem -t ${threads} ${assembly} ${pair12} ${pair22} | samtools sort -@ 4 -o ./coverage/${prefix}.bwamem.sorted.bam -
bwa mem -t ${threads} ${assembly} ${pair12} ${pair22} | samtools sort -@ 4 -o ./coverage/${prefix}.bwamem.sorted.bam

## get the coverage
bedtools genomecov -d -split -ibam ./coverage/${prefix}.bwamem.sorted.bam | gzip > ./coverage/${prefix}.bwamem.sorted.cov.tsv.gz

## calculate the median across the whole genome using the exact basepair
medianSR=$( zcat ./coverage/${prefix}.bwamem.sorted.cov.tsv.gz | awk '{if($3 != "0") print $3}' | sort -n | awk '{ a[i++]=$1} END{x=int((i+1)/2); if(x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' )
##calculate the binned median coverage and normalise each bin value by the genome wide median coverage
echo "contig;start;end;coverage_abs;coverage_norm" | tr ';' '\t' > ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.bwamem.coverage_normalised.tsv
zcat ./coverage/${prefix}.bwamem.sorted.cov.tsv.gz  | awk '{print $1"\t"$2"\t"$2"\t"$3}' | bedtools sort |\
bedtools map -b - -a ./coverage/${prefix}.${window2}kbwindow_${slide2}kbslide.bed -c 4 -o median | awk -v median="$medianSR" '{print $0"\t"$4/median}' >> ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.bwamem.coverage_normalised.tsv

if [[ $cleanup == "yes" ]]
then
##remove the alignment file due to size
rm ./coverage/${prefix}.bwamem.sorted.bam
##remove the coverage files looking at everybase pair due to size
rm ./coverage/${prefix}.bwamem.sorted.cov.tsv.gz
##remove index files
rm ${assembly}.amb
rm ${assembly}.ann
rm ${assembly}.bwt
rm ${assembly}.pac
rm ${assembly}.sa
fi

##set a full path for the shortread data to be read into the R script
SRpath=$( realpath ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.bwamem.coverage_normalised.tsv )

echo "$(date +%H:%M) ########## Step 6b: Plotting coverage"
##covert paths in R script to those relevant for the coverage outputs generated here and then export the plots 
cat ${LRSRcoverageRpath} | sed "s|PATHTOSRCOVERAGE|${SRpath}|" | sed "s|PATHTOLRCOVERAGE|${LRpath}|" | sed "s|PATHTOOUTPUT|./coverage/${prefix}|" > ./coverage/${prefix}.coverage_plots.R
Rscript ./coverage/${prefix}.coverage_plots.R
else 

echo "$(date +%H:%M) ########## Step 6b: Plotting coverage"
##covert paths in R script to those relevant for the coverage outputs generated here and then export the plots 
cat ${LRcoverageRpath} | sed "s|PATHTOLRCOVERAGE|${LRpath}|" | sed "s|PATHTOOUTPUT|./coverage/${prefix}|" > ./coverage/${prefix}.coverage_plots.R
Rscript ./coverage/${prefix}.coverage_plots.R
fi



########################## TELOMERALITY ##########################
echo "$(date +%H:%M) ########## Step 7: Running Telomere search"
mkdir ./telomerality
#telomererepeat="TTAGGG"
##can also label all regions with the canonical telomeric repeat
##use seqkit to locate the position of the canonical repeat then merge all those locations with a buffer of 7bp incase one repeat is off and export a bed file
repeatsize=$( echo $telomererepeat | awk '{print length($1)+1}' )
minrepeatsize=$( echo $telomererepeat | awk '{print length($1)*2}' )
echo "contig;start;end;sense" | tr ';' '\t' > ./telomerality/telomeres.bed
seqkit locate -j ${threads} --ignore-case -p "${telomererepeat}" ${assembly} | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d ${repeatsize} -c 4 -o distinct -i - | awk -v minrepeatsize="$minrepeatsize" '{if($3-$2 > minrepeatsize) print}' | awk -F "," '{print $1}' >> ./telomerality/telomeres.bed

##now summarise per chromosome if the ends have telomeric repeats
##the classification for presence will if the identified telomeric region is at max 75% its own length away from the assembled end
##e.g. if above 100bp of telomeric repeats was identified, this region has to be within 0.75*100bp from a contig end for that end to be considered capped by a telomere
##the classificatio for absence will be if there were no telomeres identified within 5kb of an end
##edge will be beginning or end
##closest_coords will be coordinates for the region identified
##distance_to_edge will be the shortest distance of the repeat to an edge
##classified will be categorical telomeric/distant/absent
echo "contig;edge;distance_to_edge;closest_coords;classified" | tr ';' '\t' > ./telomerality/telomeres.classification.tsv

cat ${assembly}.fai | cut -f1-2 | while read line
do
contig=$( echo "${line}" | awk '{print $1}' )
size=$( echo "${line}" | awk '{print $2}' )
cat ./telomerality/telomeres.bed | awk -v contig="$contig" '{if($1 == contig) print}' | head -n1 | awk -v contig="$contig" '{if($1 ==contig && $2 <= 0.75*($3-$2)) {print contig"\tbeginning\t"$2"\t"$2"-"$3"\ttelomeric"; end="noneed"} else if($1 ==contig && $2 > 0.75*($3-$2) && $2 < 5000) {print contig"\tbeginning\t"$2"\t"$2"-"$3"\tdistant"; end="noneed"}} END{if(end != "noneed") print contig"\tbeginning\tNA\tNA\tabsent"}'
cat ./telomerality/telomeres.bed | awk -v contig="$contig" '{if($1 == contig) print}' | tail -n1 | awk -v contig="$contig" -v size="$size" '{if($1 ==contig && (size-$3) <= 0.75*($3-$2)) {print contig"\tend\t"(size-$3)"\t"$2"-"$3"\ttelomeric"; end="noneed"} else if($1 ==contig && (size-$3) > 0.75*($3-$2) && (size-$3) < 5000) {print contig"\tend\t"(size-$3)"\t"$2"-"$3"\tdistant"; end="noneed"}} END{if(end != "noneed") print contig"\tend\tNA\tNA\tabsent"}'
done >> ./telomerality/telomeres.classification.tsv


########################## LAST CLEAN-UP ##########################
## just removing some index files, links etc
if [[ $cleanup == "yes" ]]
then
rm ${assembly}.fai
fi

########################## SUMMARY STATS ##########################

echo "$(date +%H:%M) ########## Step 8: Generating summary statistics file"
### Create the header for the summary stats file
if [[ ${sequences} == "scaffolds" ]]
then
echo "prefix;assembly;quast_#scaffolds;quast_#scaffolds>10kb;quast_assembly_size;quast_assembly_N50;quast_assembly_N90;quast_largest_contig;BUSCO_db;BUSCO_total;BUSCO_complete;BUSCO_complete_single;BUSCO_fragmented;BUSCO_missing;merqury_kmer_completeness(%);merqury_qv(phred);CRAQ_R-AQI(%);CRAQ_S-AQI(%);coverage_normal(%);telomeric_ends;telomeric_ends(%);T2T_scaffolds" | tr ';' '\t' >  summary_stats.tsv
fi

if [[ ${sequences} == "contigs" ]]
then
echo "prefix;assembly;quast_#contigs;quast_#contigs>10kb;quast_assembly_size;quast_assembly_N50;quast_assembly_N90;quast_largest_contig;BUSCO_db;BUSCO_total;BUSCO_complete;BUSCO_complete_single;BUSCO_fragmented;BUSCO_missing;merqury_kmer_completeness(%);merqury_qv(phred);CRAQ_R-AQI(%);CRAQ_S-AQI(%);coverage_normal(%);telomeric_ends;telomeric_ends(%);T2T_contigs" | tr ';' '\t' >  summary_stats.tsv
fi

##assign all the stats variables


##assembly name without the suffix (can be the prefix too if not set)
assembly2=$( echo $assembly | awk -F "/" '{print $NF}' | sed 's/\.fasta\.gz$//' | sed 's/\.fa\.gz$//' | sed 's/\.fasta$//' | sed 's/\.fa$//' | sed 's/\.fna$//' )

##QUAST
## get the information of interest out of the summary files (number of contigs; contigs > 10kb; sum size; N50, N90; largest contig) and place in variable "quaststat"
quaststat=$( cat ./quast/report.tsv | awk -F "\t" '{if(NR == 5 || NR == 14 || NR == 15 || NR == 16 || NR == 18 || NR == 19) all=all";"$NF} END{print all}' | sed 's/^;//' | tr ';' '\t' | awk '{print $2"\t"$1"\t"$4"\t"$5"\t"$6"\t"$3}' )


## BUSCO
## get the information of interest out of the summary file (busco_db;total;complete;complete_singlecopy;fragmented;missing) and place in variable "buscostat"
buscostat=$( cat ./busco/short_summary.specific.*.busco.txt | grep "The lineage\|Complete BUSCOs\|Complete and single\|Fragmented BUSCOs\|Missing BUSCOs\|Total BUSCO" | sed 's/# The lineage dataset is: //g' | awk '{line=line";"$1} END{print line}' | sed 's/^;//' | tr ';' '\t' | awk '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5}' )


## Merqury
## get the information of interest out of the two summary files (the kmer completeness percentage and the phred value for error rate) and place in variable "merqurystat"
completeness=$( cat ./merqury/${prefix}.merqury.completeness.stats | awk '{print $5}' )
phredval=$( cat ./merqury/${prefix}.merqury.qv | awk '{print $4}' )
##combine both
merqurystat=$( echo "${completeness};${phredval}" | tr ';' '\t' )


## CRAQ
## get the information of interest out of the summary file (just the final percentages for local and structural concordance) and place in variable "craqstat"
craqstat=$( cat ./craq/runAQI_out/out_final.Report | head -n3 | tail -n1 | awk -F "\t" '{print $6$7}' | tr '(' '\t' | tr ')' '\t' | awk '{print $2"\t"$4}' )

## Coverage
## calculate the proportion of the windows (sits in for genome proportion) with a median coverage less than times the standard deviation away from the median
##can just use the normalised median and therefore within 2SD from 1
covSD=$( cat ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv | awk '{print $5}' | awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' )
covstat=$( cat ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv |awk -v covSD="$covSD" '{if($5 > (1+(2*covSD)) || $5 < (1-(2*covSD))) {deviation=deviation+1}} END{print (1-(deviation/NR))*100}' )
##spit out a file with all the regions considered outside the given 2SD range
echo "contig;start;end;coverage_abs;coverage_norm" | tr ';' '\t' > ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.outside_2SD.tsv
cat ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.tsv |awk -v covSD="$covSD" '{if($5 >= (1+(2*covSD)) || $5 <= (1-(2*covSD))) {print}}' >> ./coverage/${prefix}.${window2}kbwindow_${slide2}kbsliding.minimap.coverage_normalised.outside_2SD.tsv



## telomerality
##get the number of telomere ends and as a percentage of contig ends
telomericends=$( cat ./telomerality/telomeres.classification.tsv | awk '{if($5 == "telomeric") {sum=sum+1; n++} else {n++}} END{if(sum>0) {print sum"\t"(sum/n)*100} else {print "0\t0"}}' )
##get the number of contigs with telomeres at both ends
t2t=$( cat ./telomerality/telomeres.classification.tsv | cut -f1 | sort -u | while read contig; do cat ./telomerality/telomeres.classification.tsv | awk -v contig="$contig" '{if($1 == contig && $5 == "telomeric") sum=sum+1} END{if(sum==2) print}'; done | wc -l )
##combine the two
telomeralitystat=$(  echo "${telomericends};${t2t}" | tr ';' '\t' )


##spit out all the stats and save them to the summary stats file
echo "${prefix};${assembly2};${quaststat};${buscostat};${merqurystat};${craqstat};${covstat};${telomeralitystat}" | tr ';' '\t' >> summary_stats.tsv

echo "$(date +%H:%M) ######### Summary stats can be found here ${output}/summary_stats.tsv"
echo "################ All complete; thanks for using PAQman"

