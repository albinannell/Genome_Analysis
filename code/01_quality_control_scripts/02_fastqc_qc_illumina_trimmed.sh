#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -J fastqc_qc
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/01_preprocessing/02_qc_illumina_trimmed/fastqc.%j.out

# Load modules
module load bioinfo-tools
module load FastQC

# Define variables
INPUT_FILES="/home/alan4480/genome_analysis/Genome_Analysis/analyses/01_preprocessing/02_illumina_trimmed/*.fastq.gz"  
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/01_preprocessing/02_qc_illumina_trimmed"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run FastQC on all input files
fastqc -t 4 -o $OUTPUT_DIR $INPUT_FILES

echo "Quality control completed!"

