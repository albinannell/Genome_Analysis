# Load necessary packages
suppressPackageStartupMessages({
  library(ggplot2)
})

# Set input folder and output folder (script directory)
input_dir <- "/home/alan4480/genome_analysis/Genome_Analysis/analyses/05_RNA/04_read_counting"
output_dir <- getwd()  # The folder where this script is run from

# List of count files
count_files <- c(
  "BH_ERR1797972_counts.txt",
  "BH_ERR1797973_counts.txt",
  "BH_ERR1797974_counts.txt",
  "Serum_ERR1797969_counts.txt",
  "Serum_ERR1797970_counts.txt",
  "Serum_ERR1797971_counts.txt"
)

# Loop over each file
for (file in count_files) {
  # Construct full path
  filepath <- file.path(input_dir, file)
  
  # Read the data (assumes tab-separated with no header)
  count_data <- read.table(filepath, header = FALSE, col.names = c("Gene", "Count"))
  
  # Remove any special HTSeq summary rows if present
  count_data <- subset(count_data, !grepl("^__", count_data$Gene))
  
  # Plot histogram
  p <- ggplot(count_data, aes(x = Count)) +
    geom_histogram(binwidth = 100, fill = "#2c7fb8", color = "black") +
    ggtitle(paste("Read Counts per Gene:", file)) +
    xlab("Read Count") +
    ylab("Number of Genes") +
    theme_minimal()
  
  # Save the plot
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(file), "_histogram.pdf"))
  ggsave(output_file, plot = p, width = 8, height = 6)
}

