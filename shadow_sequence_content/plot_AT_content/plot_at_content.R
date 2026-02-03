# =========================================
# GC Content Side-by-Side Boxplot per Species
# =========================================

library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

# -------------------------
# Paths
# -------------------------
COMBINED_CSV <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_sequence_content//plot_AT_content//boxplot_stats.csv"
OUTPUT_PLOT <- "C://Users//admin//Desktop//Coding//exon_shadow_exitron//shadow_sequence_content//plot_AT_content//AT_content_sidebyside.pdf"

# -------------------------
# Read the data
# -------------------------
df <- read_csv(COMBINED_CSV, show_col_types = FALSE)

# -------------------------
# Fixed species order and labels
# -------------------------
SPECIES_ORDER <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df$species <- factor(df$species, levels = SPECIES_ORDER, ordered = TRUE)
df$region <- factor(df$region, levels = c("Upstream", "Downstream"))

# -------------------------
# Wrap species names (~10 char per line)
# -------------------------
df$species_wrapped <- sapply(as.character(df$species), function(x) str_wrap(x, width = 10))

# -------------------------
# Numeric x-position for dodging boxes
# -------------------------
# Each species gets one number, then we dodge upstream/downstream around it
df$xpos <- as.numeric(df$species) + ifelse(df$region == "Upstream", -0.15, 0.15)

# Offset for N labels above max
n_offset <- 6

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(x = xpos, group = interaction(species, region), fill = region)) +
  
  # Boxplot using precomputed stats
  geom_boxplot(
    stat = "identity",
    aes(
      ymin = min, lower = q1, middle = median, upper = q3, ymax = max
    ),
    width = 0.25
  ) +
  
  # Red whisker caps
  geom_segment(aes(x = xpos - 0.125, xend = xpos + 0.125, y = min, yend = min),
               color = "red", linewidth = 1) +
  geom_segment(aes(x = xpos - 0.125, xend = xpos + 0.125, y = max, yend = max),
               color = "red", linewidth = 1) +
  
  # Median line (green)
  geom_segment(aes(x = xpos-0.125, xend = xpos+0.125 , y = median, yend = median),
               color = "darkgreen", linewidth = 0.1) +
  
  # Mean line (blue dotted)
  geom_segment(aes(x = xpos - 0.125, xend = xpos+0.125, y = mean, yend = mean),
               color = "blue", linetype = "dotted", linewidth = 0.5) +
  
  # Median/mean/min/max labels
  geom_text(aes(y = median, label = paste0(round(median, 2))),
            color = "darkgreen", vjust = -0.6, size = 2.5) +
  geom_text(aes(y = mean, label = paste0(round(mean, 2))),
            color = "blue", vjust = 1.2, size = 2.5) +
  geom_text(aes(y = min, label = paste0(round(min, 2))),
            color = "red", vjust = 1.5, size = 2.8) +
  geom_text(aes(y = max, label = paste0(round(max, 2))),
            color = "red", vjust = -0.6, size = 2.8) +
  
  # N values above max
  geom_text(aes(y = max + n_offset, label = paste0("N=", count)),
            color = "black", vjust = 0, size = 3, fontface="bold") +
  
  # X-axis: species names centered
  scale_x_continuous(
    breaks = 1:length(SPECIES_ORDER),
    labels = sapply(SPECIES_ORDER, function(x) str_wrap(x, width = 10))
  ) +
  
  # Fill colors for upstream/downstream
  scale_fill_manual(values = c("Upstream" = "pink", "Downstream" = "lightgray")) +
  
  labs(
    title = "AT Content Across Species (Upstream vs Downstream)",
    x = "Species",
    y = "AT Content (%)",
    fill = "Region"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )

# Dummy legend data
legend_df <- data.frame(
  statistic = c("Mean", "Median"),
  x = c(1, 2), y = c(0, 0)   # y values don't matter
)

# Add invisible geoms just to create legend
p <- p +
  geom_line(
    data = legend_df,
    aes(x = x, y = y, color = statistic, linetype = statistic), 
    inherit.aes = FALSE,
    size = 1,
    alpha = 1  # invisible on plot
  ) +
  scale_color_manual(
    name = "Statistic",
    values = c("Mean" = "blue", "Median" = "darkgreen")
  ) +
  scale_linetype_manual(
    name = "Statistic",
    values = c("Mean" = "dotted", "Median" = "solid")
  )

print(p)


# -------------------------
# Save plot
# -------------------------
ggsave(OUTPUT_PLOT, plot = p, width = 17, height = 7, dpi = 300)
