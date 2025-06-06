#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J trimmomatic_trim
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/01_preprocessing/02_illumina_trimmed/trimmomatic.%j.out

# Load modules
module load bioinfo-tools
module load trimmomatic/0.39

# Define variables
RAW_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/genomics/Illumina"
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/01_preprocessing/02_illumina_trimmed"
ADAPTERS="/sw/apps/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa"

# Input files
R1="$RAW_DIR/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz"
R2="$RAW_DIR/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz"

# Output files
mkdir -p $OUTPUT_DIR
OUT_P1="$OUTPUT_DIR/trimmed_R1_paired.fastq.gz"
OUT_UP1="$OUTPUT_DIR/trimmed_R1_unpaired.fastq.gz"
OUT_P2="$OUTPUT_DIR/trimmed_R2_paired.fastq.gz"
OUT_UP2="$OUTPUT_DIR/trimmed_R2_unpaired.fastq.gz"

# Run Trimmomatic
trimmomatic PE -threads 4 -phred33 \
  $R1 $R2 \
  $OUT_P1 $OUT_UP1 \
  $OUT_P2 $OUT_UP2 \
  ILLUMINACLIP:$ADAPTERS:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Trimming completed!"

