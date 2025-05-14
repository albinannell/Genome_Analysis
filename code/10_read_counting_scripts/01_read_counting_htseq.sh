#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J htseq_count
#SBATCH --mail-type=ALL
#SBATCH --mail-user=albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting/htseq.%j.out
#SBATCH --error=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting/htseq.%j.err

# Load necessary modules
module load bioinfo-tools
module load htseq/2.0.2

# Set paths
BAMDIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/03_mapping"
GFF="/home/alan4480/genome_analysis/Genome_Analysis/analyses/04_annotation/pacbio_annotation_cleaned.gff"
OUTDIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting"

mkdir -p $OUTDIR

# Loop through BAM files and run htseq-count
for BAM in $BAMDIR/*_sorted.bam; do
    SAMPLE=$(basename $BAM _sorted.bam)
    echo "Counting reads for $SAMPLE..."
    htseq-count \
        --format bam \
        --order pos \
        --stranded no \
        --type CDS \
        --idattr ID \
        $BAM $GFF > $OUTDIR/${SAMPLE}_counts.txt
done

