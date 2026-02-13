# -----------------------------
# KDE Plot Generator (TRUE FINAL VERSION)
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript kde_plotter.R <input_folder> <output_folder>")
}

input_folder <- args[1]
output_folder <- args[2]

# -----------------------------
# Smart Reader (handles tabs / spaces / mixed garbage)
# -----------------------------
smart_read <- function(file, header = FALSE) {

  df <- tryCatch(
    read.table(file, sep = "\t", header = header,
               fill = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if(is.null(df) || ncol(df) == 1) { # nolint
    df <- read.table(file, header = header,
                     fill = TRUE, stringsAsFactors = FALSE)
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

  # Valley between 0.2–0.3
  low_point_x <- NULL
  valley_idx <- which(dens$x >= 0.2 & dens$x <= 0.3)

  if (length(valley_idx) > 0) {
    low_idx <- valley_idx[which.min(dens$y[valley_idx])]
    low_point_x <- dens$x[low_idx]
  }

  # Second peak before valley
  second_peak_x <- NULL
  local_max_idx <- find_local_maxima(dens$y)

  if(!is.null(low_point_x)) {

    candidate_idx <- local_max_idx[dens$x[local_max_idx] < low_point_x]

    if(length(candidate_idx) > 0) {
      second_peak_x <- dens$x[candidate_idx[
        which.max(dens$y[candidate_idx])
      ]]
    }
  }

  # Peak between 0.3–0.4 (SCFR only)
  mid_peak_x <- NULL

  if(detect_mid_peak) {

    mid_idx <- local_max_idx[
      dens$x[local_max_idx] >= 0.3 &
      dens$x[local_max_idx] <= 0.4
    ]

    if(length(mid_idx) > 0) {
      mid_peak_x <- dens$x[mid_idx[
        which.max(dens$y[mid_idx])
      ]]
    }
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
    geom_vline(xintercept = peaks$peak_x,
               color = "red", linetype = "dashed")

  if(!is.null(peaks$second_peak_x))
    p <- p + geom_vline(xintercept = peaks$second_peak_x,
                        color = "blue", linetype = "dashed")

  if(!is.null(peaks$low_point_x))
    p <- p + geom_vline(xintercept = peaks$low_point_x,
                        color = "green", linetype = "dashed")

  if(!is.null(peaks$mid_peak_x))
    p <- p + geom_vline(xintercept = peaks$mid_peak_x,
                        color = "orange", linetype = "dashed")

  # X-axis ticks
  x_breaks <- unique(c(pretty(values),
                       peaks$peak_x,
                       peaks$second_peak_x,
                       peaks$low_point_x,
                       peaks$mid_peak_x))

  x_labels <- ifelse(
    round(x_breaks, 2) == round(peaks$peak_x, 2),
    round(x_breaks, 2),
    round(x_breaks, 3)
  )

  p +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggtitle(paste(species, "-", label, "Distribution")) +
    xlab("Frequency") +
    ylab("Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# -----------------------------
# File scanning
# -----------------------------
peak_files <- list.files(
  input_folder,
  pattern = "_with_peaks\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

scfrs_files <- list.files(
  file.path(input_folder, "scfrs"),
  pattern = "\\.tsv$",
  full.names = TRUE
)

scfrs_map <- setNames(scfrs_files,
                      sub("\\.tsv$", "", basename(scfrs_files)))

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
  # MAIN Frequencies Plot
  # -----------------------------
  df <- smart_read(file, header = TRUE)

  if(!"Frequencies" %in% colnames(df))
    stop("Missing Frequencies column in: ", file)

  freq_values <- as.numeric(unlist(
    strsplit(as.character(df$Frequencies), ";")
  ))

  main_plot <- build_plot(freq_values, species, region)

  main_pdf <- file.path(
    species_folder,
    paste0(species, "_", region, "_hist.pdf")
  )

  ggsave(main_pdf, plot = main_plot, width = 8, height = 6)

  cat("Saved:", main_pdf, "\n")

  # -----------------------------
  # SCFR Plot (SEPARATE)
  # -----------------------------
  if(species %in% names(scfrs_map)) {

    scfrs_file <- scfrs_map[[species]]
    scfrs_df <- smart_read(scfrs_file)

    cat("SCFR columns detected:", scfrs_file,
        "->", ncol(scfrs_df), "\n")

    if(ncol(scfrs_df) >= 7) {

      scfr_values <- as.numeric(unlist(
        strsplit(as.character(scfrs_df[[7]]), ";")
      ))

      scfr_plot <- build_plot(
        scfr_values,
        species,
        "SCFR",
        detect_mid_peak = TRUE
      )

      scfr_pdf <- file.path(
        species_folder,
        paste0(species, "_SCFR_hist.pdf")
      )

      ggsave(scfr_pdf, plot = scfr_plot, width = 8, height = 6)

      cat("Saved:", scfr_pdf, "\n")

    } else {
      cat("SCFR unusable:", scfrs_file, "\n")
    }
  }
}