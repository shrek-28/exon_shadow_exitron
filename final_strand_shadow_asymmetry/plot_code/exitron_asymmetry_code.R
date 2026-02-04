library(tidyr)
library(dplyr)
library(ggplot2)

df <- read.csv("C://Users//admin//Desktop//Coding//exon_shadow_exitron//final_strand_shadow_asymmetry//intermed_data//exitron_asymmetry_summary.csv", check.names = FALSE)

# -----------------------------
# 2. Long format
# -----------------------------
df_long <- df %>%
  pivot_longer(
    cols = -dataset,
    names_to = "species",
    values_to = "asymmetry"
  )

# -----------------------------
# 3. LOCK SPECIES ORDER (CRITICAL)
# -----------------------------
species_order <- c(
  "Gibbon",
  "Gorilla",
  "Human",
  "Bonobo",
  "Chimpanzee",
  "Bornean Orangutan",
  "Sumatran Orangutan"
)

df_long$species <- factor(
  df_long$species,
  levels = species_order
)

# -----------------------------
# 4. Apply line breaks AFTER factor locking
# -----------------------------
df_long$species_label <- factor(
  ifelse(
    grepl("Orangutan", df_long$species),
    str_replace(df_long$species, " ", "\n"),
    as.character(df_long$species)
  ),
  levels = c(
    "Gibbon",
    "Gorilla",
    "Human",
    "Bonobo",
    "Chimpanzee",
    "Bornean\nOrangutan",
    "Sumatran\nOrangutan"
  )
)

# -----------------------------
# 5. Plot
# -----------------------------
ggplot(df_long, aes(x = species_label, y = asymmetry, group = dataset, color = dataset)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  # STRICT y-axis
  scale_y_continuous(limits = c(-1, 1)) +
  
  # Labels
  labs(
    title = "Exitron Asymmetry Across Primate Species",
    x = "Species",
    y = "Asymmetry",
    color = "Dataset"
  ) +
  
  # Theme customization (EDIT FREELY)
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 0,
      face = "bold",
      size = 11, color="black"
    ),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(
      size = 13,      # ← X-axis label size
      face = "bold"   # optional
    ),
    axis.title.y = element_text(
      size = 13,      # ← Y-axis label size
      face = "bold"   # optional
    ),
    legend.title = element_blank()
  )

ggsave("C://Users//admin//Desktop//Coding//exon_shadow_exitron//final_strand_shadow_asymmetry//plots//exitron_asymmetry_plot.pdf", width=10, dpi=300, height=7)