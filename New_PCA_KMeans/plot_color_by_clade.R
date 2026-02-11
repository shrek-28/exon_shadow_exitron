#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(rlang)
})

# -----------------------------------
# ARGUMENTS
# -----------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_pca_biplot.R <input_folder> <output_folder>")
}

input_folder  <- args[1]
output_folder <- args[2]

scores_file    <- file.path(input_folder, "pca_scores.tsv")
loadings_file  <- file.path(input_folder, "pca_loadings.tsv")
variance_file  <- file.path(input_folder, "explained_variance.tsv")
cluster_file   <- file.path(output_folder, "cluster_assignments.tsv")

# -----------------------------------
# CHECK FILES
# -----------------------------------

for (f in c(scores_file, loadings_file, variance_file, cluster_file)) {
  if (!file.exists(f)) stop(paste("Missing file:", f))
}

# -----------------------------------
# LOAD DATA
# -----------------------------------

scores   <- read_tsv(scores_file, show_col_types = FALSE)
loadings <- read_tsv(loadings_file, show_col_types = FALSE)
variance <- read_tsv(variance_file, show_col_types = FALSE)
clusters <- read_tsv(cluster_file, show_col_types = FALSE)

# -----------------------------------
# MERGE + CLADE
# -----------------------------------

plot_df <- scores %>%
  inner_join(clusters, by = "sequence") %>%
  mutate(Clade = sub(".*::([^:]+):.*", "\\1", sequence))

# -----------------------------------
# PCs
# -----------------------------------

pc_cols  <- paste0("PC", 1:5)
pc_pairs <- combn(pc_cols, 2, simplify = FALSE)
arrow_scale <- 5
label_offset <- 1.25

# -----------------------------------
# LOOP OVER PC PAIRS - SEPARATE PDFs
# -----------------------------------

for (pair in pc_pairs) {

  x_pc <- pair[1]
  y_pc <- pair[2]

  # ---- Variance explained for axis labels ----
  var_x <- variance$variance_explained[variance$PC == x_pc]
  var_y <- variance$variance_explained[variance$PC == y_pc]

  # ---- Top 5 loadings ----
  top_loadings <- loadings %>%
    select(codon, all_of(x_pc), all_of(y_pc)) %>%
    mutate(magnitude = abs(.data[[x_pc]]) + abs(.data[[y_pc]])) %>%
    arrange(desc(magnitude)) %>%
    slice(1:5)

  # ---- PDF FILE ----
  pdf_file <- file.path(output_folder,
                        paste0("pca_clade_", x_pc, "_", y_pc, ".pdf"))

  pdf(pdf_file, width = 12, height = 6)

  # ---- PLOT ----
  p <- ggplot(plot_df, aes(x = !!sym(x_pc), y = !!sym(y_pc), color = Clade)) +
    geom_point(size = 2, alpha = 0.8) +

    geom_segment(
      data = top_loadings,
      aes(x = 0, y = 0,
          xend = !!sym(x_pc) * arrow_scale,
          yend = !!sym(y_pc) * arrow_scale),
      arrow = arrow(length = unit(0.2, "cm")),
      inherit.aes = FALSE
    ) +

    geom_text(
      data = top_loadings,
      aes(x = !!sym(x_pc) * (arrow_scale + label_offset),
          y = !!sym(y_pc) * (arrow_scale + label_offset),
          label = !!sym("codon")),
      inherit.aes = FALSE,
      size = 1.5
    ) +

    labs(
      title = paste(x_pc, "vs", y_pc),
      x = paste0(x_pc, " (", round(var_x, 2), "%)"),
      y = paste0(y_pc, " (", round(var_y, 2), "%)"),
      color = "Clade"
    ) +

    theme_minimal() +
    theme(
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

  print(p)
  dev.off()
}

cat("PCA biplots with top 5 codon loadings and variance annotated saved as separate PDFs (pca_clade_...).\n")
