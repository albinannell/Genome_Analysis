#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library("DESeq2")
  library("ggplot2")
})

# Define input and output directories
input_dir <- "/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting"
output_dir <- "/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/05_diff_expression"

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define sample files and conditions
sample_files <- c(
  "BH_ERR1797972_counts.txt",
  "BH_ERR1797973_counts.txt",
  "BH_ERR1797974_counts.txt",
  "Serum_ERR1797969_counts.txt",
  "Serum_ERR1797970_counts.txt",
  "Serum_ERR1797971_counts.txt"
)

sample_names <- c("BH_1", "BH_2", "BH_3", "Serum_1", "Serum_2", "Serum_3")
conditions <- c("BH", "BH", "BH", "Serum", "Serum", "Serum")

# Build sample table
sample_table <- data.frame(
  sampleName = sample_names,
  fileName = sample_files,
  condition = factor(conditions, levels = c("BH", "Serum"))
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_table,
  directory = input_dir,
  design = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results: Serum vs. BH
res <- results(dds, contrast = c("condition", "Serum", "BH"))
res_ordered <- res[order(res$padj), ]

# Save results
# Save cleaned results as CSV (selected columns only)
res_table <- as.data.frame(res_ordered)[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(res_table,
          file = file.path(output_dir, "deseq2_results.csv"),
          row.names = TRUE)  # `row.names = TRUE` to keep gene names

# MA plot
pdf(file.path(output_dir, "MA_plot.pdf"))
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))
dev.off()

# PCA plot with equal axis limits
rld <- rlog(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +  # Equal axis scaling
  theme_minimal()

ggsave(file.path(output_dir, "PCA_plot_equal_axes.pdf"), pca_plot)

# Volcano plot
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)

# Basic filtering: remove NA padj
volcano_data <- volcano_data[!is.na(volcano_data$padj), ]

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot") +
  theme_minimal()

ggsave(file.path(output_dir, "volcano_plot.pdf"), volcano_plot)


