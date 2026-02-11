#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(cluster)
  library(fpc)
})

# -------------------------
# ARGUMENTS
# -------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript kmeans_diagnostic.R <pca_output_folder>")
}

input_folder <- args[1]
pca_file <- file.path(input_folder, "pca_scores.tsv")

if (!file.exists(pca_file)) {
  stop("pca_scores.tsv not found in input folder.")
}

# -------------------------
# LOAD DATA
# -------------------------

pca_scores <- read_tsv(pca_file, show_col_types = FALSE)

if (!"sequence" %in% colnames(pca_scores)) {
  stop("Column 'sequence' not found.")
}

pc_cols <- paste0("PC", 1:10)

if (!all(pc_cols %in% colnames(pca_scores))) {
  stop("PC1–PC10 not found in file.")
}

X <- as.matrix(pca_scores[, pc_cols])

cat("====================================\n")
cat(" BASIC STRUCTURE\n")
cat("====================================\n")

cat("Number of sequences:", nrow(X), "\n")
cat("Number of PCs used:", ncol(X), "\n\n")

# -------------------------
# DUPLICATE CHECK
# -------------------------

cat("====================================\n")
cat(" DUPLICATE CHECK\n")
cat("====================================\n")

unique_rows <- nrow(unique(X))
dup_count <- nrow(X) - unique_rows

cat("Unique PCA rows:", unique_rows, "\n")
cat("Number of duplicated PCA rows:", dup_count, "\n\n")

# -------------------------
# VARIANCE CHECK
# -------------------------

cat("====================================\n")
cat(" PC VARIANCE\n")
cat("====================================\n")

pc_variances <- apply(X, 2, var)
print(pc_variances)
cat("\n")

# -------------------------
# NA / INFINITE CHECK
# -------------------------

cat("====================================\n")
cat(" DATA QUALITY CHECK\n")
cat("====================================\n")

cat("Any NA values:", any(is.na(X)), "\n")
cat("Any Infinite values:", any(is.infinite(X)), "\n\n")

# -------------------------
# TEST KMEANS (k=3)
# -------------------------

cat("====================================\n")
cat(" TEST KMEANS (k=3)\n")
cat("====================================\n")

set.seed(42)
km_test <- kmeans(X, centers = 3, nstart = 50)

cat("Cluster sizes:\n")
print(table(km_test$cluster))
cat("\n")

# -------------------------
# SILHOUETTE + DBI
# -------------------------

cat("====================================\n")
cat(" SILHOUETTE & DBI (k=3)\n")
cat("====================================\n")

dmat <- dist(X)

sil_test <- mean(silhouette(km_test$cluster, dmat)[,3])
cat("Silhouette score:", sil_test, "\n")

cs <- tryCatch(
  cluster.stats(dmat, km_test$cluster),
  error = function(e) NULL
)

if (!is.null(cs) && length(cs$db) == 1) {
  cat("Davies-Bouldin Index:", cs$db, "\n")
} else {
  cat("Davies-Bouldin Index: FAILED\n")
}

cat("\nDiagnostics complete.\n")