# -----------------------------
# Violin Plot Generator (All regions + SCFR)
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript violin_plotter.R <input_folder> <output_folder>")
}

input_folder  <- args[1]
output_folder <- args[2]

# -----------------------------
# Smart reader (tabs or spaces)
# -----------------------------
smart_read <- function(file, header=FALSE){
  df <- tryCatch(read.table(file, sep="\t", header=header, fill=TRUE, stringsAsFactors=FALSE),
                 error=function(e) NULL)
  if(is.null(df) || ncol(df)==1){
    df <- read.table(file, header=header, fill=TRUE, stringsAsFactors=FALSE)
  }
  return(df)
}

# -----------------------------
# Build violin plot
# -----------------------------
build_violin <- function(values, species, label){
  values <- values[!is.na(values)]
  if(length(values)==0) return(NULL)
  
  log_values <- log(values)
  
  original_median <- median(values)
  original_mean   <- mean(values)
  original_min    <- min(values)
  original_max    <- max(values)
  
  df <- data.frame(value=log_values)
  
  p <- ggplot(df, aes(x="", y = value)) +
    geom_violin(fill="lightblue", alpha=0.5) +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    geom_hline(yintercept=log(original_median+1), color="red", linetype="dashed") +
    geom_hline(yintercept=log(original_mean+1), color="green", linetype="dashed") +
    geom_hline(yintercept=log(original_min+1), color="purple", linetype="dashed") +
    geom_hline(yintercept=log(original_max+1), color="orange", linetype="dashed") +
    annotate("text", x=1.1, y=log(original_median+1), label=paste("median:",original_median), color="red", hjust=0) +
    annotate("text", x=1.1, y=log(original_mean+1), label=paste("mean:",round(original_mean,2)), color="green", hjust=0) +
    annotate("text", x=1.1, y=log(original_min+1), label=paste("min:",original_min), color="purple", hjust=0) +
    annotate("text", x=1.1, y=log(original_max+1), label=paste("max:",original_max), color="orange", hjust=0) +
    ggtitle(paste(species, "-", label, "Violin Plot")) +
    ylab("Log(Value)") + xlab("") +
    theme_minimal()
  
  return(p)
}

# -----------------------------
# Helper to get exact file
# -----------------------------
get_file <- function(folder, pattern){
  files <- list.files(folder, pattern=pattern, full.names=TRUE)
  if(length(files)==0) return(NULL)
  return(files[1])
}

# -----------------------------
# Scan species
# -----------------------------
all_species <- unique(c(
  sub("upstream_with_peaks\\.tsv$","",list.files(file.path(input_folder,"upstream"),pattern="_upstream_with_peaks\\.tsv$")),
  sub("downstream_with_peaks\\.tsv$","",list.files(file.path(input_folder,"downstream"),pattern="_downstream_with_peaks\\.tsv$")),
  sub("_exitron_with_peaks\\.tsv$","",list.files(file.path(input_folder,"exitron"),pattern="_exitron_with_peaks\\.tsv$")),
  sub("\\.tsv$","",list.files(file.path(input_folder,"scfrs"),pattern="\\.tsv$"))
))

for(species in all_species){
  upstream_file   <- get_file(file.path(input_folder,"upstream"), paste0("^",species,"upstream_with_peaks\\.tsv$"))
  downstream_file <- get_file(file.path(input_folder,"downstream"), paste0("^",species,"downstream_with_peaks\\.tsv$"))
  exitron_file    <- get_file(file.path(input_folder,"exitron"), paste0("^",species,"_exitron_with_peaks\\.tsv$"))
  scfr_file       <- get_file(file.path(input_folder,"scfrs"), paste0("^",species,"\\.tsv$"))
  
  if(all(is.null(upstream_file), is.null(downstream_file), is.null(exitron_file), is.null(scfr_file))) next
  
  # create one folder per species
  species_folder <- file.path(output_folder, species)
  if(!dir.exists(species_folder)) dir.create(species_folder, recursive=TRUE)
  
  # ---------------- Violin plots for upstream/downstream/exitron
  peak_files <- list(
    upstream=list(file=upstream_file,label="upstream"),
    downstream=list(file=downstream_file,label="downstream"),
    exitron=list(file=exitron_file,label="exitron")
  )
  
  for(pf in peak_files){
    if(!is.null(pf$file)){
      df <- smart_read(pf$file, header=TRUE)
      if("Frequencies" %in% colnames(df)){
        values <- as.numeric(unlist(strsplit(as.character(df$Frequencies), ";")))
        plot <- build_violin(values, species, pf$label)
        if(!is.null(plot)){
          ggsave(file.path(species_folder, paste0(species,"_",pf$label,"_violin.pdf")), plot=plot, width=6, height=6)
          cat("Saved violin plot:", species, pf$label,"\n")
        }
      }
    }
  }
  
  # ---------------- SCFR violin plot
  if(!is.null(scfr_file)){
    scfr_df <- smart_read(scfr_file, header=FALSE)
    if(ncol(scfr_df)>=1){
      values <- as.numeric(unlist(strsplit(as.character(scfr_df[[ncol(scfr_df)]]), ";")))
      plot <- build_violin(values, species, "SCFR")
      if(!is.null(plot)){
        ggsave(file.path(species_folder, paste0(species,"_SCFR_violin.pdf")), plot=plot, width=6, height=6)
        cat("Saved SCFR violin plot:", species,"\n")
      }
    }
  }
}