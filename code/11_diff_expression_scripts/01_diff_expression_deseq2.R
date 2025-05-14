# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

# Load DESeq2 library
library("DESeq2")

# List all *_counts.txt files in input directory
count_files <- list.files(input_dir, pattern = "_counts\\.txt$", full.names = TRUE)

# Extract sample names from filenames
sample_names <- basename(count_files)

# Create sample table
samples <- data.frame(
  sampleName = sample_names,
  fileName = count_files,
  condition = ifelse(grepl("BH", sample_names), "BH", "Serum")
)
samples$condition <- factor(samples$condition)

# Build DESeq2 dataset from HTSeq-count files
dds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(resOrdered), file = file.path(output_dir, "DESeq2_results.csv"))

# Save MA plot
pdf(file.path(output_dir, "ma_plot.pdf"))
plotMA(res, main = "DESeq2 MA-plot", ylim = c(-5, 5))
dev.off()

