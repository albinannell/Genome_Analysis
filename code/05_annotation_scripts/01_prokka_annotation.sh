#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J prokka_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/04_annotation/prokka_annotation.%j.out

# Load modules
module load bioinfo-tools
module load prokka

# Define paths
INPUT_FASTA="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/01_pacbio_assembly/assembly.contigs.fasta"
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/04_annotation"

# Create output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Run Prokka
prokka --outdir $OUTPUT_DIR --prefix pacbio_annotation --cpus 1 $INPUT_FASTA --force $INPUT_FASTA

echo "Annotation completed!"

