#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 07:00:00
#SBATCH -J pacbio_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load canu

# Define variables
GENOME_SIZE="3m" 
INPUT_FILES="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/genomics/PacBio/*"  
OUTPUT_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/pacbio_assembly"

# Run Canu
canu -p assembly -d $OUTPUT_DIR \
     genomeSize=$GENOME_SIZE \
     -pacbio $INPUT_FILES \
     useGrid=false \
     maxThreads=4

echo "Genome assembly completed!"
