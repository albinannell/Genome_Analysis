#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J trimmomatic_trim_serum
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/02_trimming/02_trimming_RNA-Seq_Serum/trimmomatic.%j.out

# Load modules
module load bioinfo-tools
module load trimmomatic/0.39

# Define variables
RAW_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/transcriptomics"
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/02_trimming/02_trimming_RNA-Seq_Serum"
ADAPTERS="/sw/apps/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa"

# Input files
INPUT_FILES="$RAW_DIR/RNA-Seq_Serum/trim_paired*"

# Output files
mkdir -p $OUTPUT_DIR

# Run Trimmomatic for all input files
for R1 in $INPUT_FILES; do
    R2="${R1/_1/_2}"  # assuming paired-end reads, change if necessary
    OUT_P1="$OUTPUT_DIR/$(basename $R1 .gz)_trimmed_R1_paired.fastq.gz"
    OUT_UP1="$OUTPUT_DIR/$(basename $R1 .gz)_trimmed_R1_unpaired.fastq.gz"
    OUT_P2="$OUTPUT_DIR/$(basename $R2 .gz)_trimmed_R2_paired.fastq.gz"
    OUT_UP2="$OUTPUT_DIR/$(basename $R2 .gz)_trimmed_R2_unpaired.fastq.gz"

    # Run Trimmomatic
    trimmomatic PE -threads 2 -phred33 \
        $R1 $R2 \
        $OUT_P1 $OUT_UP1 \
        $OUT_P2 $OUT_UP2 \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Trimming completed!"
