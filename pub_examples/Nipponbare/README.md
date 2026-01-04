# Pipeline used to test the quality of different assembly versions of the Nipponbare rice cultivar

Primarily using data from the T2T-Nipponbare project https://doi.org/10.1016/j.molp.2023.08.003


## 0. Setup

    ##set a project name variable and create a directory for all the raw and analysed results
    project="Nipponbare"
    mkdir ${project}
    cd ${project}

    threads="16"

## 1. Download assemblies
Going to use a total of 5 assemblies for Nipponbare; 1 from the T2T project and the rest from earlier versions from NCBI including the default reference <br/>

    ##create conda environment for the download of assemblies, and both prefetch/download and compression of reads
    #conda create -n ncbi_datasets htslib conda-forge::ncbi-datasets-cli seqkit bioconda::sra-tools bioconda::longreadsum
    conda activate ncbi_datasets
    
    mkdir assemblies

### 1.A Nipponbare assemblies

These assemblies contain all semi-contiguous releases of Nipponbare assemblies; including the canonical ref (GCA_001433935.1) and the gold standard T2T (GCA_034140825.1 ) <br/>
    
    ##download all from NCBI
    datasets download genome accession --assembly-source GenBank GCA_003865235.1 GCA_051403585.1 GCA_000005425.2 GCA_001433935.1 GCA_034140825.1
    unzip ncbi_dataset.zip
    rm ncbi_dataset.zip
    ls ncbi_dataset/data/ | grep -v json | while read genome
    do
    genome2=$( echo $genome | sed 's/_//g' | awk -F "." '{print $1}'  )
    cat ncbi_dataset/data/$genome/$genome*.fna | bgzip --threads ${threads} > assemblies/${genome2}.fa
    done
    mv ncbi_dataset/data/*jso* ./
    rm -r ncbi_dataset
    rm README.md
    rm md5sum.txt

## 2. Download a set of reads
Downloaded both the PacBio HiFi and ONT reads from the T2T assembly project

    mkdir reads

### 2.A PacBio HiFi reads
    
    ##pacbio download (using prefetch due to size and number of files)
    mkdir reads/pacbio
    
    prefetch --max-size 100G SRR25241090
    fasterq-dump -e ${threads} SRR25241090/SRR25241090.sra
    rm -r SRR25241090
    
    ##compress all output together
    mv SRR25241090.fastq reads/pacbio/SRR25241090.pacbio.fq
    bgzip --threads ${threads} reads/pacbio/SRR25241090.pacbio.fq

    ##get stats on the dataset quickly
    longreadsum fq -t ${threads} -i reads/pacbio/SRR25241090.pacbio.fq.gz -o reads/pacbio/SRR25241090.stats

    ##longreadsum output
    #

    ##so we have the expected ~XXX coverage and XXkb read length    

### 2.B ONT reads

    ##pacbio download (using prefetch due to size and number of files)
    mkdir reads/ont
    
    prefetch --max-size 100G SRR25241091
    fasterq-dump -e ${threads} SRR25241091/SRR25241091.sra
    rm -r SRR25241091
    
    ##compress all output together
    mv SRR25241091.fastq reads/ont/SRR25241091.ont.fq
    bgzip --threads ${threads} reads/ont/SRR25241091.ont.fq

    ##get stats on the dataset quickly
    longreadsum fq -t ${threads} -i reads/ont/SRR25241091.ont.fq.gz -o reads/ont/SRR25241091.stats

    ##longreadsum output
    #

    ##so we have the expected ~XXX coverage and XXkb read length  
    
    ##get out of environment for getting the raw reads
    conda deactivate
    

## 3. Run PAQman on all assemblies
The uses PAQman (see github READme for installation/usage instructions)

    #conda create -n paqman samtobam::paqman
    conda activate paqman
    
    ##PAQman options specific for humans is only the busco database ('-b tetrapoda') (didn't use a more specific dataset as the number of busco genes grows substatially)
    busco="poaceae"
    ##telomeric repeat is the default for humans
    ##loop through each assembly and run PAQman
    ls assemblies/*.fa.gz | while read assembly
    do
    sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
    time paqman.sh -a ${assembly} -l reads/pacbio/SRR25241090.pacbio.fq.gz -t ${threads} -b ${busco} -p ${sample} -o ${project}/${sample}_paqman
    done


## 4. Run PAQplots on the four resulting summary stats files
The uses PAQman (see github READme for installation/usage instructions)

    ##create a txt file with a list of paths to the summary files
    ls *_paqman/summary_stats.tsv > list_of_summary_files.txt
    ##run paqplots for the comparisons
    time paqplots.sh -l list_of_summary_files.txt -p ${project} -o ${project}_paqplot
  
