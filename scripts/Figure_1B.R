# Figure 1B: Number of Sequences per Coverage Category per Sample

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
})

# Load filtered data
filtered_data <- readRDS("pub_scripts/data/filtered_data_wwlong.rds")

# Calculate coverage and categorize into bins
plot_data <- filtered_data %>%
  mutate(
    coverage = covered_bases / reference_length,
    cov_category = case_when(
      coverage >= 0.9 ~ "Over 90%",
      coverage >= 0.5 & coverage < 0.9 ~ "50% to 90%",
      coverage >= 0.1 & coverage < 0.5 ~ "10% to 50%",
      TRUE ~ "Less"
    )
  ) %>%
  filter(cov_category %in% c("10% to 50%", "50% to 90%", "Over 90%")) %>%
  mutate(
    cov_category = factor(
      cov_category, 
      levels = c("10% to 50%", "50% to 90%", "Over 90%")
    )
  ) %>%
  group_by(sample_ID, cov_category) %>%
  summarise(seqs = n_distinct(sequence_name), .groups = "drop")

# Calculate median sequences per coverage category
median_category <- plot_data %>%
  group_by(cov_category) %>%
  summarise(median_seqs = median(seqs, na.rm = TRUE), .groups = "drop")

# Create plot
custom_scale <- viridis(3, option = "D")

p <- ggplot(plot_data, aes(x = cov_category, y = seqs, fill = cov_category)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "grey20") +
  geom_jitter(aes(color = cov_category), alpha = 0.4, width = 0.04, size = 1.5) +
  scale_fill_manual(values = custom_scale) +
  scale_color_manual(values = custom_scale) +
  geom_text(
    data = median_category, 
    aes(x = cov_category, y = median_seqs, label = median_seqs), 
    vjust = -0.5, 
    size = 3,
    color = "black"
  ) +
  ylim(0, max(plot_data$seqs, na.rm = TRUE) + 10) +
  labs(
    x = "Coverage Category",
    y = "# Sequences per Sample",
    title = "Number of Sequences per Coverage Category per Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks.length = unit(5, "pt"),
    axis.ticks = element_line(color = "black")
  )

print(p)
