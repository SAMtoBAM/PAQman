# Pipeline used to test the quality of 2 assembly versions of the CHM13 human cell line


## 0. Setup

    ##set a project name variable and create a directory for all the raw and analysed results
    project="human"
    mkdir ${project}

## 1. Download assemblies
This download used both `wget` <br/>
The list of assemblies was manually determined using the most complete and current (v2) and an earlier version prior to all the important additions to close many gaps (v0.7)

    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
    mv chm13v2.0.fa.gz ${project}/CHM13v2.fa.gz

    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v0.7.fasta.gz
    mv chm13.draft_v0.7.fasta.gz ${project}/CHM13v0_7.fa.gz
 
## 2. Download a set of raw Oxford nanopore reads
This download uses `wget` <br/>
The dataset was manually determined from the detailed website _https://github.com/marbl/CHM13_ (trying 

    wget https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/rel3.fastq.gz
    ##DOWNSAMPLING????? XXXXXXXXXXXXXXXXXX

## 3. Run PAQman on all assemblies
The uses PAQman (see github READme for installation/usage instructions)

    ##PAQman options specific for humans is only the busco database ('-b tetrapoda') (didn't use a more specific dataset as the number of busco genes grows substatially)
    busco="tetrapoda"
    ##ability for this repeat to detect telomeric ends was manually validated
    ##loop through each assembly and run PAQman
    ls ${project}/*.fa.gz | while read assembly
    do
    sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
    time paqman.sh -a ${project}/${sample}.fa.gz -l ${project}/XXXXXXXXXXXXXXXXXX.fq.gz -t 16 -b ${busco} -p ${sample} -o ${project}/${sample}_paqman
    done


## 4. Run PAQplots on the four resulting summary stats files
The uses PAQman (see github READme for installation/usage instructions)

    ##create a txt file with a list of paths to the summary files
    ls ${project}/*_paqman/summary_file.tsv > ${project}/list_of_summary_files.txt
    ##run paqplots for the comparisons
    time paqplots.sh -l ${project}/list_of_summary_files.txt -p ${project} -o ${project}/${project}_paqplot
  
