# Pipeline used to test the assembly quality of 5 long-read assemblies of the reference strain S288c of _Saccharomyces cerevisiae_ using PAQman


## 0. Setup

    ##set a project name variable and create a directory for all the raw and analysed results
    project="S288c"
    ##create project directory for all data
    mkdir ${project}
    cd ${project}

    ##set variable for 16 threads
    threads="16"

    ##set up and activate conda environment
    ##inlcudes paqman plus tools/mibs for the download of assemblies, and both prefetch/download and compression of reads
    conda create -n paqman htslib conda-forge::ncbi-datasets-cli seqkit bioconda::sra-tools bioconda::longreadsum conda-forge::pigz samtobam::paqman
    conda activate paqman

## 1. Download assemblies
The S288c assemblies were manually picked (GCA_000146045, GCA_002057635, GCA_016858165, GCA_022626425, GCA_902192305) <br/>
These 5 assemblies were all assembled with long read technology

    mkdir assemblies/
    
    ##download all data
    datasets download genome accession GCA_000146045 GCA_002057635 GCA_016858165 GCA_022626425 GCA_902192305
    ##unzip the data
    unzip ncbi_dataset.zip
    ##remove the zipped form
    rm ncbi_dataset.zip
    ##rename the assemblies as just the ncbi GCA assession
    ##and move into the project folder (and compress)
    ls ncbi_dataset/data/ | grep -v json | while read genome
    do
    genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
    cat ncbi_dataset/data/$genome/$genome*.fna | pigz -p ${threads} > ${project}/assemblies/$genome2.fa.gz
    done
 
## 2. Download a set of raw Oxford nanopore reads
The dataset was manually determined, selecting a recent, high coverage and reasonably long read dataset (SRR17374240) <br/>
_Note: This read dataset was used to assemble GCA_022626425; which is used in this evaluation_

    mkdir reads
    
    ##set the SRR data to download as a variable
    SRR="SRR17374240"
    ##download the fastq file
    fasterq-dump ${SRR}

    ##downsample the excessive ~800X coverage to 100X
    ##use 12Mb as estimated genome size and therefore 100X this
    #filtlong -t 1200000000 --length_weight 5 ${SRR}.fastq | pigz -p ${threads} > reads/${SRR}.filtlong100x.fq.gz
    rasusa reads -b 1200000000 ${SRR}.fastq | pigz -p ${threads} > reads/${SRR}.rasusa100x.fq.gz

    ##remove the full set of reads just taking up space
    rm ${SRR}.fastq


## 3. Run PAQman on all assemblies
The uses PAQman (see github READme for installation/usage instructions)

    ##PAQman options specific for cerevisiae include the busco database ('-b saccharomycetaceae') and the telomeric repeat (-r GGTGTG)
    busco="saccharomycetaceae"
    repeat="GGTGTG"
    ##ability for this repeat to detect telomeric ends was manually validated
    ##loop through each assembly and run PAQman
    ls assemblies/*.fa.gz | while read assembly
    do
    sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
    time paqman.sh -a assemblies/${sample}.fa.gz -l reads/${SRR}.rasusa100x.fq.gz -t ${threads} -b ${busco} -p ${sample} -o ${sample}_paqman -r ${repeat}
    done


## 4. Run PAQplots on the four resulting summary stats files
The uses PAQman (see github READme for installation/usage instructions)

    ##create a txt file with a list of paths to the summary files
    ls *_paqman/summary_stats.tsv > list_of_summary_files.txt
    ##run paqplots for the comparisons
    paqplots.sh -l list_of_summary_files.txt -p ${project} -o ${project}_paqplot


 <p align="center" >
    <img src="https://github.com/SAMtoBAM/PAQman/blob/main/figures/figure2.paqplot_example.S288c_alone.svg" width=70%>
</p>
