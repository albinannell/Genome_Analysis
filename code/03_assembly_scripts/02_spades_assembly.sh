#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J spades_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/02_illumina_nanopore_assembly/spades.%j.out

# Load modules
module load bioinfo-tools
module load spades

# Define directories and files
ILLUMINA_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/genomics/Illumina"
NANOPORE_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/genomics/Nanopore"
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/02_illumina_nanopore_assembly"

ILLUMINA_R1="$ILLUMINA_DIR/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz"
ILLUMINA_R2="$ILLUMINA_DIR/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz"
NANOPORE="$NANOPORE_DIR/E745_all.fasta.gz"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run SPAdes with Illumina paired-end and Nanopore reads for hybrid assembly with single k-mer size (e.g., 55)
spades.py --pe1-1 $ILLUMINA_R1 --pe1-2 $ILLUMINA_R2 \
          --nanopore $NANOPORE \
          --k-mer-size 55 \
          -o $OUTPUT_DIR \
          -t 2

echo "SPAdes genome assembly completed!"

