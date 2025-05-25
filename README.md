# Genome_Analysis
The following folder provides the pipeline for doing several bioinformatic analyses based on data from the bacteria E. faecium.

The data are from the article "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum" by Zhang et. al that used transcriptome profiling (RNA-Seq) and a transposon mutant library (Tn-Seq) to find genes that contribute to growth of E. faecium in human serum. The pipeline provided here aims to similairly find genes that contribute to the growth in E. faecium in human serum. 

In order to achieve that this folder contains data, scripts for running the analyses and results and outputs from the anaylses performed. These three categories are divided into three subfolders "data", "code" and "analyses". 

The main pipeline is as follows:
1. Genome Assembly 
2. Assembly Evaluation
3. Annotation
4. Differential Expression Analysis

In more detail the pipeline is meant to be run in the following order:
1. Genome Assembly using PacBio data and Canu
2. Assembly Evaluation using QUAST
3. Illumina DNA reads QC using FastQC
4. Illumina reads preprocessing using Trimmomatic
5. Illumina DNA reads QC using FastQC
6. Illumina/Nanopore genome assembly using Spades
7. Assembly Comparison using MUMmerplot
8. Assembly Evaluation using QUAST
9. Annotation using Prokka
10. Synteny comparison of PacBio genome assembly and reference using ACT
11. Illumina RNA reads QC using FastQC
12. Illumina reads preprocessing using Trimmomatic
13. Illumina RNA reads QC using FastQC
14. RNA mapping using BWA/SAM-tools
15. Read counting using Htseq
16. Differential expression analysis using Deseq2
