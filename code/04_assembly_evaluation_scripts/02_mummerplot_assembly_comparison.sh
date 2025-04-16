#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J mummer_comparison
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/03_assembly_evaluation/02_mummerplot_comparison/%x.%j.out

# Load modules
module load bioinfo-tools
module load MUMmer

# Define input files
REF="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/01_pacbio_assembly/assembly.contigs.fasta"
QUERY="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/02_illumina_nanopore_assembly/contigs.fasta"

# Output directory
OUTDIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/03_assembly_evaluation/02_mummerplot_comparison"
mkdir -p $OUTDIR
cd $OUTDIR

# Run nucmer
nucmer --prefix=spades_vs_canu $REF $QUERY

# Generate alignment plot with contig delimiters
mummerplot --fat --filter --postscript --layout -R $REF -Q $QUERY -p spades_vs_canu spades_vs_canu.delta

echo "Assembly comparison completed! Plot with contig delimiters saved in $OUTDIR"

