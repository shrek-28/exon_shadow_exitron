#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

# -----------------------------
# Command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2){
  stop("Usage: Rscript boxplot_metric.R <input_file> <output_pdf>")
}

input_file <- args[1]
output_pdf <- args[2]

# -----------------------------
# Load data
# -----------------------------
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure Metric is numeric
df$Metric <- as.numeric(df$Metric)

# -----------------------------
# Generate boxplot
# -----------------------------
p <- ggplot(df, aes(y = Metric)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.color = "red") +
  ylab("Metric") +
  ggtitle("Distribution of SCFR Metrics") +
  theme_minimal()

# -----------------------------
# Save to PDF
# -----------------------------
ggsave(output_pdf, plot = p, width = 6, height = 4)
cat("Boxplot saved to:", output_pdf, "\n")
