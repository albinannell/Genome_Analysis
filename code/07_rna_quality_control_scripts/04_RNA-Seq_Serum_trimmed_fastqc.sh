#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J fastqc_rna_trimmed_Serum_qc
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/01_qc/04_qc_RNA-Seq_Serum_trimmed/fastqc.%j.out

# Load modules
module load bioinfo-tools
module load FastQC

# Define variables
INPUT_FILES="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/02_trimming/02_trimming_RNA-Seq_Serum/*_paired.fastq.gz"
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/01_qc/04_qc_RNA-Seq_Serum_trimmed"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run FastQC on all input files
fastqc -t 2 -o $OUTPUT_DIR $INPUT_FILES

echo "RNA quality control completed!"
