#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(cluster)
  library(clusterSim)
})

# -----------------------------
# ARGUMENTS
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript robust_kmeans_rule_based.R <input_folder> <output_folder>")

input_folder <- args[1]
output_folder <- args[2]

pca_file <- file.path(input_folder, "pca_scores.tsv")
if (!file.exists(pca_file)) stop("pca_scores.tsv not found in input folder.")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# -----------------------------
# LOAD DATA & PREPARE MATRIX
# -----------------------------
pca_scores <- read_tsv(pca_file, show_col_types = FALSE)
if (!"sequence" %in% colnames(pca_scores)) stop("Column 'sequence' not found.")

num_pcs <- min(20, sum(grepl("^PC", colnames(pca_scores))))
pc_cols <- paste0("PC", 1:num_pcs)
X <- as.matrix(pca_scores[, pc_cols])
n_seq <- nrow(X)

# -----------------------------
# K RANGE
# -----------------------------
k_values <- 2:10
sil_values <- numeric(length(k_values))
dbi_values <- numeric(length(k_values))
wcss_values <- numeric(length(k_values))

# -----------------------------
# SAFE NEAREST-CENTROID ASSIGNMENT
# -----------------------------
assign_nearest <- function(points, centers) {
  points <- as.matrix(points)
  centers <- as.matrix(centers)
  n_points <- nrow(points)
  n_centers <- nrow(centers)
  
  dists <- matrix(0, nrow = n_points, ncol = n_centers)
  for (i in 1:n_centers) {
    dists[, i] <- rowSums((points - centers[i, ])^2)
  }
  max.col(-dists)
}

# -----------------------------
# SILHOUETTE SAMPLING PARAMETERS
# -----------------------------
sil_threshold <- 5000   # Dataset size above which silhouette is sampled
max_sample_size <- 1000 # Maximum number of points to sample

# -----------------------------
# OPTIMIZATION LOOP
# -----------------------------
set.seed(42)
for (i in seq_along(k_values)) {
  k <- k_values[i]
  
  # Mini-batch KMeans for speed
  batch_idx <- sample(1:n_seq, min(5000, n_seq))
  km <- kmeans(X[batch_idx, , drop = FALSE], centers = k, nstart = 30)
  
  # Assign all points to nearest centroid
  cluster_assign <- assign_nearest(X, km$centers)
  
  # Compute WCSS and DBI
  wcss_values[i] <- sum(rowSums((X - km$centers[cluster_assign, , drop = FALSE])^2))
  dbi_values[i] <- clusterSim::index.DB(X, cluster_assign)$DB
  
  # -----------------------------
  # SILHOUETTE (DYNAMIC SAMPLING)
  # -----------------------------
  if (n_seq > sil_threshold) {
    sample_size <- min(max_sample_size, n_seq)
    sil_idx <- sample(1:n_seq, sample_size)
  } else {
    sil_idx <- 1:n_seq
  }
  
  sample_clusters <- cluster_assign[sil_idx]
  cluster_counts <- table(sample_clusters)
  valid_idx <- sil_idx[sample_clusters %in% names(cluster_counts[cluster_counts > 1])]
  
  if (length(valid_idx) > 0) {
    sil_sample <- silhouette(cluster_assign[valid_idx], dist(X[valid_idx, , drop = FALSE]))
    sil_values[i] <- mean(sil_sample[, 3])
  } else {
    sil_values[i] <- NA
  }
}

# -----------------------------
# ELBOW DETECTION (CURVATURE)
# -----------------------------
curvature <- rep(NA, length(wcss_values))
for (i in 2:(length(wcss_values) - 1)) {
  curvature[i] <- wcss_values[i - 1] - 2 * wcss_values[i] + wcss_values[i + 1]
}

# -----------------------------
# RULE-BASED HYBRID K SELECTION
# -----------------------------
valid_idx <- which(!is.na(sil_values) & sil_values > 0)

if (length(valid_idx) > 0) {
  # Step 1: pick k with lowest DBI among valid silhouettes
  candidate_idx <- valid_idx[which.min(dbi_values[valid_idx])]
  
  # Step 2: if multiple tie, pick highest curvature
  if (length(candidate_idx) > 1) {
    candidate_idx <- candidate_idx[which.max(curvature[candidate_idx])]
  }
  
  best_k <- k_values[candidate_idx]
} else {
  # fallback: pick k with highest curvature
  best_k <- k_values[which.max(curvature)]
}

# -----------------------------
# SAVE METRICS
# -----------------------------
metrics_df <- data.frame(
  k = k_values,
  Silhouette = sil_values,
  DBI = dbi_values,
  WCSS = wcss_values,
  Curvature = curvature
)
write_tsv(metrics_df, file.path(output_folder, "k_optimization_scores.tsv"))

# -----------------------------
# FINAL KMEANS ON FULL DATA
# -----------------------------
set.seed(42)
batch_idx <- sample(1:n_seq, min(5000, n_seq))
init_km <- kmeans(X[batch_idx, , drop = FALSE], centers = best_k, nstart = 50)
final_clusters <- assign_nearest(X, init_km$centers)

# -----------------------------
# SAVE CLUSTER ASSIGNMENTS
# -----------------------------
cluster_assignments <- data.frame(
  sequence = pca_scores$sequence,
  Cluster = final_clusters
)
write_tsv(cluster_assignments, file.path(output_folder, "cluster_assignments.tsv"))

# -----------------------------
# WRITE CLUSTER TXT FILES
# -----------------------------
clusters <- sort(unique(final_clusters))
for (cl in clusters) {
  seq_ids <- cluster_assignments$sequence[cluster_assignments$Cluster == cl]
  writeLines(seq_ids, file.path(output_folder, paste0("cluster_", cl, "_sequences.txt")))
}

# -----------------------------
# CLUSTER SUMMARY TSV
# -----------------------------
summary_file <- file.path(output_folder, "cluster_summary.tsv")
summary_conn <- file(summary_file, "w")
writeLines("Cluster\tn_sequences", summary_conn)
for (cl in clusters) {
  n_seq_cl <- sum(final_clusters == cl)
  writeLines(paste0(cl, "\t", n_seq_cl), summary_conn)
}
close(summary_conn)

cat("Best k selected (Rule-based hybrid: silhouette + DBI + curvature):", best_k, "\n")
cat("Stable, memory-efficient, warning-free KMeans complete.\n")
