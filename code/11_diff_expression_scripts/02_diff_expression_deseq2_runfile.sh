#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH --job-name=deseq2_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/05_diff_expression/deseq2.%j.out

# Load R and required modules
module load bioinfo-tools
module load R/4.3.1
module load R_packages/4.3.1

# Run the R script
Rscript /home/alan4480/genome_analysis/Genome_Analysis/code/11_diff_expression_scripts/01_diff_expression_deseq2.R \
  /home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting \
  /home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/05_diff_expression

