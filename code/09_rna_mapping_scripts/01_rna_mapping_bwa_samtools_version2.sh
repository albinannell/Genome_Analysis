#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J bwa_rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user albin.annell.4480@student.uu.se
#SBATCH --output=/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/03_mapping/bwa_mapping.%j.out

# Load modules
module load bioinfo-tools
module load bwa
module load samtools

# Set paths
REF="/home/alan4480/genome_analysis/Genome_Analysis/analyses/02_genome_assembly/01_pacbio_assembly/assembly.contigs.fasta"
OUTDIR="/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/03_mapping"
BH_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/transcriptomics/RNA-Seq_BH"
SERUM_DIR="/home/alan4480/genome_analysis/Genome_Analysis/data/raw_data/transcriptomics/RNA-Seq_Serum"

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Index reference genome if not already indexed
if [ ! -e "${REF}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index $REF
fi

# Function to map paired reads
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

# Map BH samples
map_reads "$BH_DIR/trim_paired_ERR1797972_pass_1.fastq.gz" "$BH_DIR/trim_paired_ERR1797972_pass_2.fastq.gz" "BH_ERR1797972"
map_reads "$BH_DIR/trim_paired_ERR1797973_pass_1.fastq.gz" "$BH_DIR/trim_paired_ERR1797973_pass_2.fastq.gz" "BH_ERR1797973"
map_reads "$BH_DIR/trim_paired_ERR1797974_pass_1.fastq.gz" "$BH_DIR/trim_paired_ERR1797974_pass_2.fastq.gz" "BH_ERR1797974"

# Map Serum samples
map_reads "$SERUM_DIR/trim_paired_ERR1797969_pass_1.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797969_pass_2.fastq.gz" "Serum_ERR1797969"
map_reads "$SERUM_DIR/trim_paired_ERR1797970_pass_1.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797970_pass_2.fastq.gz" "Serum_ERR1797970"
map_reads "$SERUM_DIR/trim_paired_ERR1797971_pass_1.fastq.gz" "$SERUM_DIR/trim_paired_ERR1797971_pass_2.fastq.gz" "Serum_ERR1797971"

echo "RNA mapping completed!"

