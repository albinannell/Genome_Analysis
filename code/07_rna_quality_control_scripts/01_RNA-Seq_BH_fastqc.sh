#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J fastqc_rna_qc
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/01_qc/01_qc_RNA-Seq_BH/fastqc.%j.out

# Load modules
module load bioinfo-tools
module load FastQC

# Define variables
INPUT_FILES="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/transcriptomics/RNA-Seq_BH/trim_paired*"  
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/01_qc/01_qc_RNA-Seq_BH"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run FastQC on all input files
fastqc -t 2 -o $OUTPUT_DIR $INPUT_FILES

echo "RNA quality control completed!"

