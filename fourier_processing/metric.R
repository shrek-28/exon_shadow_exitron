# =========================================
# KDE + Row-level Metric Plot Generator
# =========================================
suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------
# Command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) stop("Usage: Rscript kde_plotter.R <input_folder> <output_folder>")

input_folder <- args[1]
output_folder <- args[2]

# -----------------------------
# Smart Reader (handles tabs / spaces)
# -----------------------------
smart_read <- function(file, header = FALSE, sep = "\t") {
  df <- tryCatch(
    read.table(file, sep = sep, header = header, fill = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if(is.null(df) || ncol(df) == 1) {
    df <- read.table(file, header = header, fill = TRUE, stringsAsFactors = FALSE)
  }
  return(df)
}

# -----------------------------
# Filename parser
# -----------------------------
parse_filename <- function(file_path) {
  fname <- basename(file_path)
  if (grepl("upstream", fname, ignore.case = TRUE)) {
    region <- "upstream"
    species <- sub("upstream_with_peaks\\.tsv$", "", fname)
  } else if(grepl("downstream", fname, ignore.case = TRUE)) {
    region <- "downstream"
    species <- sub("downstream_with_peaks\\.tsv$", "", fname)
  } else if(grepl("exitron", fname, ignore.case = TRUE)) {
    region <- "exitron"
    species <- sub("_?exitron_with_peaks\\.tsv$", "", fname)
  } else {
    stop("Unknown region in file: ", file_path)
  }
  list(species = species, region = region)
}

# -----------------------------
# Local maxima detector
# -----------------------------
find_local_maxima <- function(y) {
  which(diff(sign(diff(y))) == -2) + 1
}

# -----------------------------
# Peak computation
# -----------------------------
compute_peaks <- function(values, detect_mid_peak = FALSE) {
  dens <- density(values, na.rm = TRUE)
  peak_x <- dens$x[which.max(dens$y)]
  low_point_x <- NULL
  valley_idx <- which(dens$x >= 0.2 & dens$x <= 0.3)
  if(length(valley_idx) > 0) low_point_x <- dens$x[valley_idx[which.min(dens$y[valley_idx])]]
  second_peak_x <- NULL
  local_max_idx <- find_local_maxima(dens$y)
  if(!is.null(low_point_x)) {
    candidate_idx <- local_max_idx[dens$x[local_max_idx] < low_point_x]
    if(length(candidate_idx) > 0) second_peak_x <- dens$x[candidate_idx[which.max(dens$y[candidate_idx])]]
  }
  mid_peak_x <- NULL
  if(detect_mid_peak) {
    mid_idx <- local_max_idx[dens$x[local_max_idx] >= 0.3 & dens$x[local_max_idx] <= 0.4]
    if(length(mid_idx) > 0) mid_peak_x <- dens$x[mid_idx[which.max(dens$y[mid_idx])]]
  }
  list(
    peak_x = peak_x,
    second_peak_x = second_peak_x,
    low_point_x = low_point_x,
    mid_peak_x = mid_peak_x
  )
}

# -----------------------------
# Plot builder
# -----------------------------
build_plot <- function(values, species, label, detect_mid_peak = FALSE) {
  peaks <- compute_peaks(values, detect_mid_peak)
  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    geom_vline(xintercept = peaks$peak_x, color = "red", linetype = "dashed")
  if(!is.null(peaks$second_peak_x)) p <- p + geom_vline(xintercept = peaks$second_peak_x, color = "blue", linetype = "dashed")
  if(!is.null(peaks$low_point_x)) p <- p + geom_vline(xintercept = peaks$low_point_x, color = "green", linetype = "dashed")
  if(!is.null(peaks$mid_peak_x)) p <- p + geom_vline(xintercept = peaks$mid_peak_x, color = "orange", linetype = "dashed")
  x_breaks <- unique(c(pretty(values), peaks$peak_x, peaks$second_peak_x, peaks$low_point_x, peaks$mid_peak_x))
  x_labels <- ifelse(round(x_breaks, 2) == round(peaks$peak_x, 2), round(x_breaks, 2), round(x_breaks, 3))
  p +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggtitle(paste(species, "-", label, "Distribution")) +
    xlab("Frequency / Metric") +
    ylab("Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------
# Row-level metric computation
# -----------------------------
compute_row_metric <- function(freq_vec, mag_vec) {
  if(length(freq_vec) != length(mag_vec)) stop("Freq and Mag length mismatch")
  if(sum(mag_vec, na.rm = TRUE) == 0) return(NA)
  sum(freq_vec * mag_vec, na.rm = TRUE) / sum(mag_vec, na.rm = TRUE)
}

# -----------------------------
# File scanning
# -----------------------------
peak_files <- list.files(input_folder, pattern = "_with_peaks\\.tsv$", recursive = TRUE, full.names = TRUE)
scfrs_files <- list.files(file.path(input_folder, "scfrs"), pattern = "\\.tsv$", full.names = TRUE)
scfrs_map <- setNames(scfrs_files, sub("\\.tsv$", "", basename(scfrs_files)))

# -----------------------------
# Main loop
# -----------------------------
for(file in peak_files) {

  info <- parse_filename(file)
  species <- info$species
  region <- info$region

  species_folder <- file.path(output_folder, species)
  if(!dir.exists(species_folder)) dir.create(species_folder, recursive = TRUE)

  # -----------------------------
  # Read main frequency + magnitude file
  # -----------------------------
  df <- smart_read(file, header = TRUE)
  if(!all(c("Frequencies","Magnitudes") %in% colnames(df))) stop("Missing Frequencies or Magnitudes in: ", file)

  # Compute row-level metric
  row_metrics <- sapply(1:nrow(df), function(i) {
    freq_values <- as.numeric(unlist(strsplit(as.character(df$Frequencies[i]), ";")))
    mag_values  <- as.numeric(unlist(strsplit(as.character(df$Magnitudes[i]), ";")))
    compute_row_metric(freq_values, mag_values)
  })

  # -----------------------------
  # Plot raw frequency distribution
  # -----------------------------
  freq_values_all <- unlist(lapply(1:nrow(df), function(i) {
    as.numeric(unlist(strsplit(as.character(df$Frequencies[i]), ";")))
  }))
  main_plot <- build_plot(freq_values_all, species, region)
  main_pdf <- file.path(species_folder, paste0(species,"_",region,"_hist.pdf"))
  ggsave(main_pdf, plot = main_plot, width = 8, height = 6)
  cat("Saved:", main_pdf, "\n")

  # -----------------------------
  # Plot metric distribution
  # -----------------------------
  metric_plot <- build_plot(row_metrics, species, paste0(region,"_metric"))
  metric_pdf <- file.path(species_folder, paste0(species,"_",region,"_metric_hist.pdf"))
  ggsave(metric_pdf, plot = metric_plot, width = 8, height = 6)
  cat("Saved:", metric_pdf, "\n")

  # -----------------------------
  # SCFRs (space-separated, no header)
  # -----------------------------
  if(species %in% names(scfrs_map)) {
    scfrs_file <- scfrs_map[[species]]
    scfrs_df <- smart_read(scfrs_file, header = FALSE, sep = " ")

    if(ncol(scfrs_df) >= 8) {
      row_metrics <- sapply(1:nrow(scfrs_df), function(i) {
        freq_values <- as.numeric(unlist(strsplit(as.character(scfrs_df[[7]][i]), " ")))
        mag_values  <- as.numeric(unlist(strsplit(as.character(scfrs_df[[8]][i]), " ")))
        compute_row_metric(freq_values, mag_values)
      })

      scfr_metric_plot <- build_plot(row_metrics, species, "SCFR_metric", detect_mid_peak = TRUE)
      scfr_metric_pdf <- file.path(species_folder, paste0(species,"_SCFR_metric_hist.pdf"))
      ggsave(scfr_metric_pdf, plot = scfr_metric_plot, width = 8, height = 6)
      cat("Saved:", scfr_metric_pdf, "\n")
    } else {
      cat("SCFR unusable:", scfrs_file, "\n")
    }
  }
}