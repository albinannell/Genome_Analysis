#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:40:00
#SBATCH -J quast_pacbio
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/03_quast_evaluation/pacbio/quast.%j.out

# Load modules
module load bioinfo-tools
module load quast

# Define variables
GENOME_SIZE="3m"
ASSEMBLY="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/pacbio_assembly/assembly.contigs.fasta"  
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/03_quast_evaluation/pacbio"

#Create output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Run QUAST
quast.py $ASSEMBLY -o $OUTPUT_DIR --gene-finding

echo "QUAST evaluation completed!"

