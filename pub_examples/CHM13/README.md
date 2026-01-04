# Pipeline used to test the quality of 3 assembly versions of the CHM13 human cell line

Primarily using data from the T2T-CHM13 github https://github.com/marbl/CHM13


## 0. Setup

    ##set a project name variable and create a directory for all the raw and analysed results
    project="CHM13"
    mkdir ${project}
    cd ${project}

    threads="16"

## 1. Download assemblies
Going to use a total of 6 assemblies for CHM13; 4 from the T2T-CHM13 project and another two that were generated earlier using long-reads and are available on NCBI <br/>
The two additional NCBI assemblies are to show a much less contiguous assembly quality relative to the latest possibilities

    ##create conda environment for the download of assemblies, and both prefetch/download and compression of reads
    #conda create -n ncbi_datasets htslib conda-forge::ncbi-datasets-cli seqkit bioconda::sra-tools
    conda activate ncbi_datasets
    
    mkdir assemblies

### 1.A T2T-CHM13 assemblies

These assemblies contain all publically released T2T versions including the most complete and current (gapless + Y chromosome) (CHM13v2; T2T-CHM13v2.0; GCA_009914755.4), The first gapless version (CHM13v1.0; T2T-CHM13v1.0; GCA_009914755.2); another with some additional improvements such as corrections and polishing (CHM13v1.1; T2T-CHM13v1.1; GCA_009914755.3) and the earliest version prior to all the important additions to close many gaps (T2T-CHM13v0.7; GCA_009914755.1) <br/>
    
    ##first get the T2T-CHM13 genomes
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
    mv chm13v2.0.fa.gz assemblies/CHM13v2.fa.gz

    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
    mv chm13.draft_v1.1.fasta.gz assemblies/CHM13v1_1.fa.gz

    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz
    mv chm13.draft_v1.0.fasta.gz assemblies/CHM13v1_0.fa.gz
    
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v0.7.fasta.gz
    mv chm13.draft_v0.7.fasta.gz assemblies/CHM13v0_7.fa.gz

    

### 2.B Other publically available CHM13 assemblies using long-reads
These two assemblies GCA_002884485.1 and GCA_000983455.2 were assembled much earlier than the T2T assemblies

    ##download two other earlier reference quality, long-read assemblies for CHM13
    datasets download genome accession --assembly-source GenBank GCA_002884485.1 GCA_000983455.2
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## 2. Download a set of raw Oxford nanopore reads
This download uses `wget` <br/>
The dataset was manually determined from the detailed website _https://github.com/marbl/CHM13/blob/master/Sequencing_data.md_ 

    mkdir reads

### 2.A PacBio HiFi reads
    
    ##pacbio download (using prefetch due to size and number of files)
    ##only able to get the ~25X coverage of 10kb HiFi reads on NCBI (the 20kb are not available for an unknown reason; only the raw subreads from the github)
    ##split into 4 different SRA submissions within https://www.ncbi.nlm.nih.gov/sra/SRX5633451/
    mkdir reads/pacbio

    ##prefetch the SRA data etc
    prefetch SRR9087597
    ##download the reads sequentially in order to help remove the temporary prefetch files after each download
    fasterq-dump -e ${threads} SRR9087597/SRR9087597.sra
    ##remove the prefetch data
    rm -r SRR9087597

    prefetch SRR9087598
    fasterq-dump -e ${threads} SRR9087598/SRR9087598.sra
    rm -r SRR9087598

    prefetch SRR9087599
    fasterq-dump -e ${threads} SRR9087599/SRR9087599.sra
    rm -r SRR9087599

    prefetch SRR9087600
    fasterq-dump -e ${threads} SRR9087600/SRR9087600.sra
    rm -r SRR9087600
    
    ##compress all output together
    cat *.fastq | bgzip --threads ${threads} > reads/pacbio/SRR9087XXX.pacbio.fq.gz
    ##remove uncompressed
    rm *.fastq

    ##get stats on the dataset quickly
    longreadsum fq -t 16 -i reads/pacbio/SRR9087XXX.pacbio.fq.gz -o ont/stats

    ##longreadsum output
    #total number of reads	6899727
    #total number of bases	75627496237
    #longest read length	21034
    #N50 read length	10859
    #mean read length	10960.94
    #median read length	10729
    #GC%	41.02

    ##so we have the expected ~24X coverage and 10kb read length
    

### 2.B ONT reads

    ##ONT fastq reads
    mkdir reads/ont
    
    ##using the latest (and easily available) set of basecalled reads, therefore the highest quality; here the 'rel8'
    ##however we don't need all 100X+; therefore we will unzip as we download, and select only reads >7.5kb in length and only until we reach ~25X which is about the same as the pacbio data
    wget -qO- https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/nanopore/rel8-guppy-5.0.7/reads.fastq.gz | gunzip -c | seqkit seq -m 7500 | seqkit head -l 75000000000 | bgzip --threads ${threads} > reads/ont/rel8.25X_subset.ont.fq.gz


        ##get stats on the dataset quickly
    longreadsum fq -t 16 -i reads/ont/rel8.25X_subset.ont.fq.gz -o ont/stats

    ##longreadsum output
    #

    ##so we have the expected ~XXX coverage and ana average XXkb read length
    
    ##get out of environment for getting the assemblies and reads
    conda deactivate
    

## 3. Run PAQman on all assemblies
The uses PAQman (see github READme for installation/usage instructions)

    #conda create -n paqman samtobam::paqman
    conda activate paqman
    
    ##PAQman options specific for humans is only the busco database ('-b tetrapoda') (didn't use a more specific dataset as the number of busco genes grows substatially)
    busco="tetrapoda"
    ##telomeric repeat is the default for humans
    ##loop through each assembly and run PAQman
    ls assemblies/*.fa.gz | while read assembly
    do
    sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
    time paqman.sh -a ${assembly} -l reads/pacbio/SRR9087XXX.pacbio.fq.gz -t ${threads} -b ${busco} -p ${sample} -o ${project}/${sample}_paqman
    done


## 4. Run PAQplots on the four resulting summary stats files
The uses PAQman (see github READme for installation/usage instructions)

    ##create a txt file with a list of paths to the summary files
    ls *_paqman/summary_stats.tsv > list_of_summary_files.txt
    ##run paqplots for the comparisons
    time paqplots.sh -l list_of_summary_files.txt -p ${project} -o ${project}_paqplot
  
