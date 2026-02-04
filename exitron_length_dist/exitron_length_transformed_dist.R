library(ggplot2)
library(dplyr)

# Read data
log_data <- read.csv("C://Users//admin//Desktop//Coding//exon_shadow_exitron//exitron_length_dist//inter_data//log_intron_lengths.csv")
stats <- read.csv("C://Users//admin//Desktop//Coding//exon_shadow_exitron//exitron_length_dist//inter_data//intron_stats.csv")

SPECIES_ORDER <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

fix_species <- function(df) {
  df %>%
    mutate(
      species = case_when(
        species == "borangutan" ~ "Bornean Orangutan",
        species == "sorangutan" ~ "Sumatran Orangutan",
        TRUE ~ tools::toTitleCase(species)
      ),
      species = factor(species, levels = SPECIES_ORDER)
    )
}

df$species_wrapped <- sapply(as.character(df$species), function(x) str_wrap(x, width = 10))

log_data <- fix_species(log_data)
stats    <- fix_species(stats)

# Compute log-transformed stats for lines
stats <- stats %>%
  mutate(log_mean = log10(orig_mean ),
         log_median = log10(orig_median),
         log_min = log10(orig_min ),
         log_max = log10(orig_max))



# Start boxplot (whiskers extend to min/max, no outliers)
p <- ggplot(log_data, aes(x=species, y=log_intron_length)) +
  geom_boxplot(fill="lightgrey", color="black", coef = Inf, width=0.5) +  # whiskers = min/max
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(color="black")) +
  ylab("Log10(Exitron Length)") +
  xlab("Species") +
  ggtitle("Exitron Length Distribution")

# Overlay lines
p <- p +
  # Median (dark green solid)
  geom_segment(data=stats,
               aes(x=as.numeric(factor(species)) - 0.25,
                   xend=as.numeric(factor(species)) + 0.25,
                   y=log_median, yend=log_median),
               color="darkgreen", size=1.2) +
  # Mean (blue dotted)
  geom_segment(data=stats,
               aes(x=as.numeric(factor(species)) - 0.25,
                   xend=as.numeric(factor(species)) + 0.25,
                   y=log_mean, yend=log_mean),
               color="blue", linetype="dotted", size=0.75) +
  # Min (red)
  geom_segment(data=stats,
               aes(x=as.numeric(factor(species)) - 0.25,
                   xend=as.numeric(factor(species)) + 0.25,
                   y=log_min, yend=log_min),
               color="red", size=1.2) +
  # Max (red)
  geom_segment(data=stats,
               aes(x=as.numeric(factor(species)) - 0.25,
                   xend=as.numeric(factor(species)) + 0.25,
                   y=log_max, yend=log_max),
               color="red", size=1.2)

# Annotate original values
p <- p +
  geom_text(data=stats,
            aes(x=species, y=log_median, label=paste0(orig_median, " bp")),
            color="darkgreen", size=3, vjust=-0.6) +
  geom_text(data=stats,
            aes(x=species, y=log_mean, label=paste0(round(orig_mean,2), " bp")),
            color="blue", size=3, vjust=-1) +
  geom_text(data=stats,
            aes(x=species, y=log_min, label=paste0(orig_min, " bp")),
            color="red", size=3, vjust=1.2) +
  geom_text(data=stats,
            aes(x=species, y=log_max, label=paste0(round(orig_max/1000,2), " kb")),
            color="red", size=3, vjust=-0.6) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

p <- p +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

p <- p +
  geom_text(
    data = stats,
    aes(x = species, y = Inf, label = paste0("N=", count)),
    vjust = 1.3,
    size = 3.5,
    fontface = "bold"
  )
  
# Dummy legend data
legend_df <- data.frame(
  statistic = c("Mean", "Median"),
  x = c(1, 2), y = c(1, 1)   # y values don't matter
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


ggsave("C://Users//admin//Desktop//Coding//exon_shadow_exitron//exitron_length_dist//exitron_length_dist_plot.pdf", plot=p, width=8, height=6)