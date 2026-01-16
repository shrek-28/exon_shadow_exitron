library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)

# ---- Data preparation ----
strand_counts <- read.csv("C:\\Users\\admin\\Documents\\internship\\strand_asymmetry_distribution\\summary_table_composite_exon.csv",
               stringsAsFactors = FALSE)

strand_counts$species <- str_wrap(strand_counts$species, width = 10)
strand_counts$species <- factor(strand_counts$species, levels = unique(strand_counts$species))

# ---- Plot with wrapped labels and small black text ----
strand_plot <- ggplot(strand_counts, aes(x = species, y = count, fill = strand)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8), 
           width = 0.6) +
  
  # Add count labels above each bar
  geom_text(aes(label = count),
            position = position_dodge(width = 0.8), 
            vjust = -0.5,  
            size = 3.25,       # smaller text
            color = "black") +
  
  labs(title = "Strand Asymmetry Analysis - Composite Exon",
       x = "Species", y = "Count") +
  scale_fill_manual(values = c("+" = "darkgrey", "-" = "#f2b949"),
                    name = "Strand") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 10), # small black text
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

strand_plot 

ggsave("C:\\Users\\admin\\Documents\\internship\\strand_asymmetry_distribution\\strand_asymmetry_all_species_analysis_composite_exon.pdf", plot = strand_plot, width = 12, height = 8, dpi = 300, bg='white')

