# Figure 1A: Percentage of Aligned Reads

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
})

# Load filtered data
filtered_data <- readRDS("pub_scripts/data/filtered_data_wwlong.rds")

# Compute % aligned per sample
plot_data <- filtered_data %>%
  group_by(sample_ID, total_filtered_reads_in_sample) %>%
  summarise(target_aligned = sum(reads_aligned, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    percent_aligned = (target_aligned / total_filtered_reads_in_sample) * 100
  )

mean_perc <- plot_data %>%
  summarise(mean_percent = mean(percent_aligned, na.rm = TRUE), .groups = "drop")

# Colors
col_black  <- "#000000"
col_yellow <- "#FFF3B0"
col_red    <- "#D32B2B"

p <- ggplot(plot_data, aes(x = "", y = percent_aligned)) +
  geom_violin(fill = col_yellow, color = col_black, alpha = 0.5, trim = FALSE) +
  geom_jitter(width = 0.06, size = 2.2, color = col_red, alpha = 0.7) +
  geom_boxplot(width = 0.08, outlier.shape = NA, fill = NA, color = col_black, linewidth = 0.4) +
  labs(title = "Percentage of Aligned Reads", x = NULL, y = "Aligned reads (%)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = col_black, fill = NA, size = 0.8),
    axis.ticks = element_line(color = col_black),
    axis.text.x = element_text(color = col_black, face = "bold"),
    axis.text.y = element_text(color = col_black),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  geom_text(data = mean_perc,
            aes(x = "", y = mean_percent + 5,
                label = paste0("Mean: ", round(mean_percent, 1), "%")),
            color = col_black, size = 3.8) +
  ylim(0, NA)

print(p)

ggsave(filename = "Figure_1A_wwlong.png", plot = p, width = 6, height = 5, dpi = 300)

cat("Saved: Figure_1A_wwlong.png\n")
