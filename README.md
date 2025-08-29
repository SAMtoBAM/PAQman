 <p align="center" >
    <img src="https://github.com/SAMtoBAM/PAQman/blob/main/logo/paqman_logo_grey.svg" width=70%>
</p>

[![Zenodo DOI](https://zenodo.org/badge/952795493.svg)](https://doi.org/10.5281/zenodo.16039705)
[![Anaconda_version](https://anaconda.org/samtobam/paqman/badges/version.svg)](https://anaconda.org/samtobam/paqman)
[![Anaconda_platforms](https://anaconda.org/samtobam/paqman/badges/platforms.svg)](https://anaconda.org/samtobam/paqman)
[![Anaconda_downloads](https://anaconda.org/samtobam/paqman/badges/downloads.svg)](https://anaconda.org/samtobam/paqman)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/paqman/badges/latest_release_date.svg)](https://anaconda.org/samtobam/paqman)


PAQman combines a set of excellent tools (and PAQman scripts) for a reference-free and comprehensive evaluation of genome assemblies <br/>

PAQman evaluates 7 important features: <br/>
Contiguity ([QUAST](http://dx.doi.org/10.1093/nar/gkad406)) <br/>
Gene Content ([BUSCO](https://doi.org/10.1093/nar/gkae987)) <br/>
Completeness ([Merqury](https://doi.org/10.1186/s13059-020-02134-9))<br/>
Error rate ([Merqury](https://doi.org/10.1186/s13059-020-02134-9))<br/>
Correctness ([CRAQ](https://doi.org/10.1038/s41467-023-42336-w))<br/>
Coverage (PAQman (with use of [bwa](https://doi.org/10.48550/arXiv.1303.3997)/[minimap2](https://doi.org/10.1093/bioinformatics/bty191) + [samtools](https://doi.org/10.1093/bioinformatics/btp352) + [bedtools](https://doi.org/10.1093/bioinformatics/btq033)))<br/>
Telomerality* (PAQman (with use of [seqtk](https://github.com/lh3/seqtk) + [bedtools](https://doi.org/10.1093/bioinformatics/btq033))) <br/>
<i>PAQman is clearly built upon the back of the tools in brackets so please cite them</i>

***

### Apptainer usage
```
docker pull ghcr.io/samtobam/paqman:latest
```

### Conda installation
```
conda config --append channels pwwang
conda install samtobam::paqman
```

### Quick run
```
paqman.sh -a path/to/assembly.fa -l path/to/long-reads.fq.gz
```


```
paqman.sh -g assembly.fa -l long-reads.fq.gz -x ont -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz -b eukaryota -w 30000 -s 10000 -r TTAGGG -p assembly -o paqman_output -c yes

Required inputs:
-a | --assembly     Genome assemly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz)
-l | --longreads    Long reads used for assembly in fastq or fasta format  (*.fa / *.fasta / *.fna / *.fastq / *.fq) and can be gzipped (*.gz)

Recommended inputs:
-x | --platform     Long-read technology to determine mapping mapping parameters. Choose between 'ont' or 'pacbio-hifi' or 'pacbio-clr' (default: ont)
-b | --buscodb      Name of BUSCO database to be used (default: eukaryota)
-t | --threads      Number of threads for tools that accept this option (default: 1)
-r | --repeat       Telomeric repeat pattern (default: TTAGGG)

Optional parameters:
-1 | --pair1        Paired end illumina reads in fastq format; first pair. Used by Merqury, CRAQ and coverage analysis (Recommended). Can be gzipped (*.gz)
-2 | --pair2        Paired end illumina reads in fastq format; second pair. Used by Merqury, CRAQ and coverage analysis (Recommended). Can be gzipped (*.gz)
-w | --window       Number of basepairs for window averaging for coverage (default: 30000)
-s | --slide        Number of basepairs for the window to slide for coverage (default: 10000)
-p | --prefix       Prefix for output (default: name of assembly file (-a) before the fasta suffix)
-o | --output       Name of output folder for all results (default: paqman_output)
-seq | --sequences	Whether or not to use scaffolds or contigs; provide 'scaffolds' to not break the assembly at N's (default: contigs)
-c | --cleanup      Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
-h | --help         Print this help message
```

***

## The Pipeline

 <p align="center" >
    <img src="https://github.com/SAMtoBAM/PAQman/blob/main/figures/paqman_schematic.svg" width=100%>
</p>


## The summary output metrics:
Although PAQman looks at 7 features, within these features are many important metrics for complete assembly evaluation <br/>
PAQman extracts some of the most informative/important and places them into a summary file <br/>
All metrics are detailed below

### Summary stats: 'summary_stats.tsv':
|Column | Header | Description |
|:---:|:---:|--------------|
| 01 | <b>prefix</b> | prefix given to the output files (-p)
| 02 | <b>assembly</b> | prefix from the fasta file without suffix (-a)
| 03 | <b>quast_#contigs</b> | Number of total contigs
| 04 | <b>quast_#contigs>10kb</b> | Number of contigs >10 kb
| 05 | <b>quast_assembly_size</b> | Total number of basepairs in assembly
| 06 | <b>quast_assembly_N50</b> | Assembly N50
| 07 | <b>quast_assembly_N90</b> | Assembly N90
| 08 | <b>quast_largest_contig</b> | Largest contig in the assembly
| 09 | <b>BUSCO_db</b> | The BUSCO database used to evaluate the assembly
| 10 | <b>BUSCO_total</b> | Total number of BUSCOs in the database
| 11 | <b>BUSCO_complete</b> | BUSCOs identified as complete
| 12 | <b>BUSCO_complete_single</b> | BUSCOs identified as complete and as a single copy
| 13 | <b>BUSCO_fragmented</b> | BUSCOs identified as fragmented
| 14 | <b>BUSCO_missing</b> | BUSCOs not identified
| 15 | <b>merqury_kmer_completeness(%)</b> | A k-mer estimation of completeness
| 16 | <b>merqury_qv(phred)</b> | A k-mer estimation of the genome wide error rate
| 17 | <b>CRAQ_R-AQI(%)</b> | Quality measure from 0-100 based on small regional errors
| 18 | <b>CRAQ_S-AQI(%)</b> | Quality measure from 0-100 based on large structural errors
| 19 | <b>coverage_normal(%)</b> | Percentage of the genome within 2*stdev of the genome wide median
| 20 | <b>telomeric_ends</b> | Number of contig ends with telomeric repeats
| 21 | <b>telomeric_ends(%)</b> | Percentage of contig ends with telomeric repeats
| 22 | <b>t2t_contigs</b> | Number of contigs with telomeric repeats at both ends


## *BUSCO dataset
To find the best appropriate BUSCO dataset for your assembly you can refer to [here](https://busco.ezlab.org/busco_userguide.html#obtain-busco) <br/>
Or run `busco --list-datasets` 


## *Telomerality:
<i>*I am using the term to describe stats about how many of the assembled contig have reached telomere sequences giving confidence of structural completeness at contig ends in repeat regions </i><br/>
Telomerality is calculated using a few simple steps specific to PAQman <br/>

1. All exact single repeats are identified (seqkit locate --ignore-case -p ${telomererepeat}) <br/>
2. Coordinates of single repeats are merged if withint 7 bp (allowing for 1 repeat to deviate mildly) <br/>
3. Only keep regions where at least two consecutive repeats were found (i.e. only keep region > 2\*repeat length) <br/>
Can find all the coordinates for telomeric regions (including interstitial) in the bed file with explanatory header: 'telomerality/telomeres.bed' <br/>
4. Contig ends are labelled in 3 ways <br/>
   <b>telomeric</b>: Coordinates for a telomeric repeat are at least within 0.75\*length from the end (e.g. a 100 bp long telomeric repeat region with within 75bp of a contig end) <br/>
   <b>distant</b>: >0.75\*length bp away but within 5kb <br/>
   <b>absent</b>: >5kb from the end or no repeats identified in contig <br/>
Can find these classifications (and coordinates/distance from edge etc) for each contig end in the tsv file with explanatory header: 'telomerality/telomeres.classification.tsv' <br/>

Note: For the option -r (--repeat); although some repeats are not exact this can still work as the detection scheme allows for inexact repeats. For example 'GGTGTG' works very well for <i>S. cerevisiae</i>, which usually is represented as T(G)*1-3.


## Comparing PAQman output across multiple assemblies

PAQman also has a tool (paqplots) to compare and analyse summary files from multiple assemblies <br/>
This makes it easier to benchmark tools and parameters using all the variables analysed  <br/>
Simply provide paqplots with a combined summary file (with the same header) or a list of paths to the summary files

```
paqplots.sh -s summary_file.tsv -p prefix -o paqplot_output
OR
paqplots.sh -l list_of_summary_files.txt -p prefix -o paqplot_output
	
Required inputs:
-s | --summary     A PAQman summary file with multiple assemblies combined and the same header
-l | --list		   A list of paths to multiple summary files

Optional parameters:
-p | --prefix       Prefix for output (default: paqplot)
-o | --output       Name of output folder for all results (default: paqplot_output)
-h | --help         Print this help message	
```

paqplot will provide two images (alongside the R scripts used to generate them) <br/>
First: Figures containing all of the raw variables compared <br/>
A radar plot (thanks ggradar!) for stats of percentages and a lollipop plot for all others (due to the wildly different scales)

<p align="center" >
    <img src="https://github.com/SAMtoBAM/PAQman/blob/main/figures/example.raw_values.svg" width=70%>
</p>

Second: Another radar plot (thanks ggradar again!) containing a subset of stats and their relative values (i.e. all stats divided by the maximum value for that stat) <br/>
In this example, all stats should be maximised except for contig count hence the PACman like shape in blue below

<p align="center" >
    <img src="https://github.com/SAMtoBAM/PAQman/blob/main/figures/example.relative_values.svg" width=50%>
</p>










