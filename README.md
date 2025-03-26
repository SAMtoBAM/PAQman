
 # Post-Assembly Quality manager ![](https://github.com/SAMtoBAM/PAQman/blob/main/paqman_logo_grey.svg)

PAQman is combines a set of tools (and a few in-house scripts) designed to comprehensively evaluate a de-novo genome assembly using a range of important metrics and without the need for a reference <br/>
<i/>Note: PAQman has been primarily designed for long-read reference quality assemblies </i>

PAQman evaluate 6 important metrics: <br/>
Contiguity; Gene Content; Completeness; Error rate; Correctness; Coverage; Telomerality*.

***

### Conda installation
```
conda install paqman
```

### Quick run
```
paqman -g genome.fa -l long-reads.fq.gz
```


```
paqman -g genome.fa -l long-reads.fq.gz -x ont -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz -b eukaryota -w 30000 -s 10000 -r TTAGGG -p genome -o paqman_output -c yes

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
-c | --cleanup      Remove a large number of files produced by each of the tools that can take up a lot of space. Choose between 'yes' or 'no' (default: 'yes')
-h | --help         Print this help message
```

***

## The output:

### Summary stats: 'summary_stats.tsv':
|Column &nbsp;&nbsp;&nbsp;| Header | Description |
|------------|--------|--------------|
| Column 01 | <b>strain</b> | prefix given to the output files (-p)
| Column 02 | <b>assembly</b> | prefix from the fasta file without suffix (-g)
| Column 03 | <b>quast_#contigs</b> | Number of total contigs
| Column 04 | <b>quast_#contigs>10kb</b> | Number of contigs >10 kb
| Column 05 | <b>quast_assembly_N50</b> | Assembly N50
| Column 06 | <b>quast_assembly_N90</b> | Assembly N90
| Column 07 | <b>quast_largest_contig</b> | Largest contig in the assembly
| Column 08 | <b>BUSCO_db</b> | The BUSCO database used to evaluate the assembly
| Column 09 | <b>BUSCO_total</b> | Total number of BUSCOs in the database
| Column 10 | <b>BUSCO_complete_single</b> | BUSCOs identified as complete and as a single copy
| Column 11 | <b>BUSCO_fragmented</b> | BUSCOs identified as fragmented
| Column 12 | <b>BUSCO_missing</b> | BUSCOs not identified
| Column 13 | <b>merqury_completeness(%)</b> | A k-mer estimation of compmeteness
| Column 14 | <b>merqury_qv(phred)</b> | A k-mer estimation of the genome wide error rate
| Column 15 | <b>CRAQ_average_CRE(%)</b> | An estimation of the total assembly without small regional errors
| Column 16 | <b>CRAQ_average_CSE(%)</b> | An estimation of the total assembly without large structural errors
| Column 17 | <b>telomeric_ends</b> | Number of contig ends with telomeric repeats
| Column 18 | <b>telomeric_ends(%)</b> | Percentage of contig ends with telomeric repeats
| Column 19 | <b>t2t_contigs</b> | Number of contigs with telomeric repeats at both ends



## *Telomerality:
<i>*I am using the term to describe stats about how many of the assembled contig have reached telomere sequences giving confidence of structural completeness at contig ends in repeat regions </i><br/>
Telomerality is calculated using a few simple steps specific to PAQman <br/>

1. All exact single repeats are identified (seqkit locate --ignore-case -p ${telomererepeat}) <br/>
2. Coordinates of single repeats are merged if withint 7 bp (allowing for 1 repeat to deviate mildly) <br/>
3. Only keep regions where at least two consecutive repeats were found (i.e. only keep region > 2\*repeat length) <br/>
Can find all the coordinates for telomeric regions (including interstitial) in the bed file with explanatory header: 'telomerality/telomeres.bed' <br/>
4. Contig ends are labelled in 3 ways <br/>
       A. 'telomeric': Coordinates for a telomeric repeat are at least within 0.75\*length from the end (e.g. a 100 bp long telomeric repeat region with within 75bp of a contig end) <br/>
       B. 'distant': >0.75\*length bp away but within 5kb <br/>
       C. 'absent': >5kb from the end or no repeats identified in contig <br/>
Can find these classifications (and coordinates/distance from edge etc) for each contig end in the tsv file with explanatory header: 'telomerality/telomeres.classification.tsv' <br/>










