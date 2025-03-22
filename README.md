# GenomEval
GenomEval is a compilation of tools designed to comprehensively evaluate a genome assembly using a range of different and important metrics


The ease of assembling a reference-quality genome has improved dramatically in the past few years <br/>
However, good practices on how to evaluate and compare the resulting assembly has lagged behind <br/>
Although there exists tools for the job, they are often used sporadically <br/>
Because of this the different and important metrics that they calculate are often overlooked making it difficult to fully evaluate an assembly and have a complete image of quality  <br/>
This is where GenomEval steps in: 

GenomeEval primarily combines a set of tools that each evaluate their own important metrics and provides you with a summary of the compiled results <br/>
Now if you ever run an assembler, GenomEval can be your next step for comprehensive genome evaluation statistics <br/>
Additionally GenomEval only requires what most people already has access to for the assembly process so there is not need to gather more data <br/>

Note: GenomEval has been designed for long-read reference quality assemblies

The 6 important metrics and the tools chosen to evaluate them:

Contiguity: Quast (great easy tool for quick asembly evaluation) <br/>
Gene content: BUSCO (measure of how well coding sequences have been assembled) <br/>
Completeness and Error rate: Merqury (k-mer based approach for both stats) <br/>
Correctness: Filtlong + CRAQ (long-read downsampling + read alignment detection of potential assembly error) <br/>
Coverage: bwa-mem/minimap2 + samtools/bedtools + genomeval-specific (alignment + processing of coverage + visualising relative coverage for each contig) <br/>
Telomerality: seqkit + genomeval-specific (identifying telomere sequences for indications of telomerality) <br/>

Note: I am using telomerality as a term to describe stats about how many of the assembled contig have reached telomere sequences giving confidence of structural completeness at contig ends in repeat regions.

```
genomeval -g genome.fa -l long-reads.fq.gz -x ont -1 illumina.R1.fq.gz -2 illumina.R2.fq.gz -b eukaryota -w 30000 -s 10000 -r TTAGGG -p genome -o genomeval_output

Required inputs:
-g -genome        Genome assemly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz)
-l --longreads    Long reads used for assembly in fastq or fasta format  (*.fa / *.fasta / *.fna / *.fastq / *.fq) and can be gzipped (*.gz)

Optional parameters:
-x --platform     Long-read technology. Choose between 'ont' or 'pacbio' (default: ont)
-1                Paired end illumina reads; first pair. Used by Merqury, CRAQ and coverage analysis
-2                Paired end illumina reads; second pair. Used by Merqury, CRAQ and coverage analysis
-b --busco        Name of BUSCO database to be used (default: eukaryota)
-w --window       Number of basepairs for window averaging for coverage (default: 30000)
-s --slide        Number of basepairs for the window to slide for coverage (default: 10000)
-r --repeat       Telomeric repeat pattern (Default: TTAGGG)
-p --prefix       Prefix for output (default: name -g before the fasta prefix)
-p --prefix       Prefix for output (default: genomeeval_output)
-t --threads      Number of threads for tools that accept this option (Default: 1)
```






