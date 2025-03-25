# GenomEval
GenomEval is combines a set of tools (and a few in-house scripts) designed to comprehensively evaluate a de-novo genome assembly using a range of important metrics and without the need for a reference <br/>
<i/>Note: GenomEval has been primarily designed for long-read reference quality assemblies </i>

The ease of assembling a reference-quality genome has improved dramatically in the past few years <br/>
However, good practices on how to evaluate and compare the resulting assembly has lagged behind <br/>
Although there exists tools for the job, they are often used sporadically <br/>
Because of this the important metrics that they calculate are often overlooked making it difficult to fully evaluate an assembly and have a complete image of quality  <br/>
This is where GenomEval steps in: 

GenomeEval combines a set of tools (and a few in-house scripts) that each evaluate their own important metrics; and provides you with a summary of the compiled results <br/>
Now if you ever run an assembler, GenomEval can be your next step for comprehensive genome evaluation statistics <br/>
The summary statistics can be used to comprehensively evaluate different assemblies for the same dataset/sample <br/>
And no reference genome is required for any of the statistics <br/>
Laslty, GenomEval only requires what most people already have access to for the assembly process, long- and (or not) short reads, so there is not need to gather more data <br/>


The 6 important metrics and the tools chosen to evaluate them:

Contiguity: Quast (great easy tool for quick asembly evaluation) <br/>
Gene content: BUSCO (uses commonly conserved genes to measure of how well coding sequences have been assembled) <br/>
Completeness and Error rate: Merqury (a k-mer based approach to determine how much of the entire genomic material has been assembled and the rate of error) <br/>
Correctness: Filtlong + CRAQ (long-read downsampling + read alignment detection of structural assembly error) <br/>
Coverage: bwa-mem/minimap2 + samtools/bedtools + genomeval-specific (alignment + processing of coverage + visualising relative coverage for each contig) <br/>
Telomerality*: seqkit + genomeval-specific (identifying telomere sequences + determining if contigs have terminal telomeric repeats) <br/>

<i>*I am using 'Telomerality' as a term to describe stats about how many of the assembled contig have reached telomere sequences giving confidence of structural completeness at contig ends in repeat regions </i>


```
genomeval -g genome.fa -l long-reads.fq.gz -x ont -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz -b eukaryota -w 30000 -s 10000 -r TTAGGG -p genome -o genomeval_output -c yes

Required inputs:
-a --assembly        Genome assemly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz)
-l --longreads    Long reads used for assembly in fastq or fasta format  (*.fa / *.fasta / *.fna / *.fastq / *.fq) and can be gzipped (*.gz)

Optional parameters:
-x --platform     Long-read technology to determine mapping mapping parameters. Choose between 'ont' or 'pacbio-hifi' or 'pacbio-clr' (default: ont)
-1 --pair1        Paired end illumina reads; first pair. Used by Merqury, CRAQ and coverage analysis (Recommended)
-2 --pair2        Paired end illumina reads; second pair. Used by Merqury, CRAQ and coverage analysis (Recommended)
-b --buscodb      Name of BUSCO database to be used (default: eukaryota)
-w --window       Number of basepairs for window averaging for coverage (default: 30000)
-s --slide        Number of basepairs for the window to slide for coverage (default: 10000)
-r --repeat       Telomeric repeat pattern (default: TTAGGG)
-p --prefix       Prefix for output (default: name of assembly file (-a) before the fasta suffix)
-o --output       Name of output folder for all results (default: genomeeval_output)
-t --threads      Number of threads for tools that accept this option (default: 1)
-c --cleanup      Remove a large number of files produced by each of the tools. Choose between 'yes' or 'no' (default: 'yes')
-h --help         Print this help message
```


## Installation and quick start

### Conda installation
```
```

### Quick run
```
```


## The output:

### Summary stats: 'summary_stats.tsv':
Column 01: 'strain': prefix given to the output files so you can easily compare across samples <br/>
Column 02: 'assembly':  prefix from the fasta file used as input so you can easily compare across assemblies <br/>
Column 03: 'quast_#contigs': Number of total contigs <br/>
Column 04: 'quast_#contigs>10kb': Number of contigs >10kb in size <br/>
Column 05: 'quast_assembly_N50': Assembly N50 <br/>
Column 06: 'quast_assembly_N90': Assembly N90 <br/>
Column 07: 'quast_largest_contig': Largest contig in the assembly <br/>
Column 08: 'BUSCO_db': The BUSCO database used to evaluate the assembly (easy to be sure of comparing BUSCO values from the same database) <br/>
Column 09: 'BUSCO_total': Total number of BUSCOs in the database <br/>
Column 10: 'BUSCO_complete_single': BUSCOs identified as complete and as a single copy <br/>
Column 11: 'BUSCO_fragmented': BUSCOs identified as fragmented <br/>
Column 12: 'BUSCO_missing': BUSCOs not identified <br/>
Column 13: 'merqury_completeness(%)': A k-mer estimation of the amount of total genomic material assembled <br/>
Column 14: 'merqury_qv(phred)': A k-mer estimation of the genome wide error rate <br/>
Column 15: 'CRAQ_average_CRE(%)': An estimation of the total assembly without any small regional errors <br/>
Column 16: 'CRAQ_average_CSE(%)': An estimation of the total assembly without any large structural errors <br/>
Column 17: 'telomeric_ends': Number of contig ends identified with telomeric repeats <br/>
Column 18: 'telomeric_ends(%)': Percentage of contig ends with telomeric repeats <br/>
Column 19: 't2t_contigs': Number of contigs with telomeric repeats at both ends <br/>



## Telomerality steps:
1. All exact single repeats are identified (seqkit locate --ignore-case -p ${telomererepeat}) <br/>
2. Coordinates of single repeats are merged if withint 7 bp (allowing for 1 repeat to deviate mildly) <br/>
3. Only keep regions where at least two consecutive repeats were found (i.e. only keep region > 2\*repeat length) <br/>
Can find all the coordinates for telomeric regions (including interstitial) in the bed file with explanatory header: 'telomerality/telomeres.bed' <br/>
4. Contig ends are labelled in 3 ways <br/>
       A. 'telomeric': Coordinates for a telomeric repeat are at least within 0.75\*length from the end (e.g. a 100 bp long telomeric repeat region with within 75bp of a contig end) <br/>
       B. 'distant': >0.75\*length bp away but within 5kb <br/>
       C. 'absent': >5kb from the end or no repeats identified in contig <br/>
Can find these classifications (and coordinates/distance from edge etc) for each contig end in the tsv file with explanatory header: 'telomerality/telomeres.classification.tsv' <br/>












