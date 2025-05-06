#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J bwa_rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/03_mapping/bwa_mapping.%j.out

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools

# Variables
REF="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/01_pacbio_assembly/assembly.contigs.fasta"
OUTDIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/03_mapping"
BH_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/02_trimming/01_trimming_RNA-Seq_BH"
SERUM_DIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/02_trimming/02_trimming_RNA-Seq_Serum"

# Create output directory
mkdir -p $OUTDIR

# Index reference genome if not already done
bwa index $REF

# Function to map paired reads with direct piping to sorted BAM
map_reads () {
    local R1=$1
    local R2=$2
    local SAMPLE=$3

    echo "Mapping $SAMPLE..."
    bwa mem -t 2 $REF $R1 $R2 | \
        samtools view -@ 2 -b - | \
        samtools sort -@ 2 -o $OUTDIR/${SAMPLE}_sorted.bam
    samtools index $OUTDIR/${SAMPLE}_sorted.bam
}

# BH samples
map_reads "$BH_DIR/trim_paired_ERR1797972_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$BH_DIR/trim_paired_ERR1797972_pass_2.fastq_trimmed_R2_paired.fastq.gz" "BH_ERR1797972"
map_reads "$BH_DIR/trim_paired_ERR1797973_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$BH_DIR/trim_paired_ERR1797973_pass_2.fastq_trimmed_R2_paired.fastq.gz" "BH_ERR1797973"
map_reads "$BH_DIR/trim_paired_ERR1797974_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$BH_DIR/trim_paired_ERR1797974_pass_2.fastq_trimmed_R2_paired.fastq.gz" "BH_ERR1797974"

# Serum samples
map_reads "$SERUM_DIR/trim_paired_ERR1797969_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797969_pass_2.fastq_trimmed_R2_paired.fastq.gz" "Serum_ERR1797969"
map_reads "$SERUM_DIR/trim_paired_ERR1797970_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797970_pass_2.fastq_trimmed_R2_paired.fastq.gz" "Serum_ERR1797970"
map_reads "$SERUM_DIR/trim_paired_ERR1797971_pass_1.fastq_trimmed_R1_paired.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797971_pass_2.fastq_trimmed_R2_paired.fastq.gz" "Serum_ERR1797971"

echo "All RNA mapping steps completed!"

