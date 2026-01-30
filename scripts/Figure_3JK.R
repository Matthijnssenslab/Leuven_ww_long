# Norovirus Strain Composition and Read Count Analysis
# Uses 500bp coverage threshold filtered data

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(tibble)
library(purrr)
library(cowplot)


# Load pre-filtered data (500bp threshold)
filtered_data <- readRDS("data/filtered_data_wwlong_500bp.rds")

all_samples <- filtered_data %>% 
  distinct(sample_ID, Date)

norwalk_data <- filtered_data %>% 
  filter(species == "Norwalk virus")

# Calculate strain-level RPKMF

strain_rpkmf <- norwalk_data %>%
  group_by(sample_ID, Date, strain) %>%
  summarise(total_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

no_detection_samples <- all_samples %>%
  anti_join(strain_rpkmf, by = c("sample_ID", "Date")) %>%
  mutate(strain = "No Norwalk Detected", total_RPKMF = 0)

combined_rpkmf <- bind_rows(strain_rpkmf, no_detection_samples)

# Identify top strains and group others
num_top_strains <- 7

top_strains <- combined_rpkmf %>%
  filter(strain != "No Norwalk Detected") %>%
  group_by(strain) %>%
  summarise(total_RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_RPKMF)) %>%
  slice_head(n = num_top_strains) %>%
  pull(strain)

combined_rpkmf <- combined_rpkmf %>%
  mutate(strain_grouped = case_when(
    strain %in% top_strains ~ strain,
    strain != "No Norwalk Detected" ~ "Other",
    TRUE ~ strain
  )) %>%
  group_by(sample_ID, Date, strain_grouped) %>%
  summarise(RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop")

# Calculate relative abundance per sample
combined_rpkmf <- combined_rpkmf %>%
  group_by(sample_ID, Date) %>%
  mutate(total_sample_RPKMF = sum(RPKMF, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_abundance = ifelse(
    total_sample_RPKMF > 0,
    (RPKMF / total_sample_RPKMF) * 100,
    ifelse(strain_grouped == "No Norwalk Detected", 100, 0)
  ))

# Create sample labels and ordering
combined_rpkmf <- combined_rpkmf %>%
  arrange(Date, sample_ID) %>%
  mutate(sample_label = paste0(sample_ID, "\n", format(Date, "%Y-%m-%d")))

sample_label_levels <- unique(combined_rpkmf$sample_label)
combined_rpkmf$sample_label <- factor(combined_rpkmf$sample_label, 
                                      levels = sample_label_levels)

# Calculate read counts per sample (for panel 1)
reads_per_sample <- filtered_data %>%
  distinct(sample_ID, Date) %>%
  left_join(
    norwalk_data %>%
      group_by(sample_ID, Date) %>%
      summarise(read_count = sum(reads_aligned, na.rm = TRUE), .groups = "drop"),
    by = c("sample_ID", "Date")
  ) %>%
  replace_na(list(read_count = 0)) %>%
  arrange(Date, sample_ID) %>%
  mutate(
    sample_label = paste0(sample_ID, "\n", format(Date, "%Y-%m-%d")),
    sample_label = factor(sample_label, levels = sample_label_levels),
    log10_reads = ifelse(read_count > 0, log10(read_count), 0)
  )

# Create x-axis tick marks (April and October of each year)
date_range <- range(combined_rpkmf$Date)
all_years <- seq(year(date_range[1]), year(date_range[2]), 1)
tick_dates <- as.Date(c(outer(all_years, c("04-01", "10-01"), paste, sep = "-")))
tick_dates <- tick_dates[tick_dates >= date_range[1] & tick_dates <= date_range[2]]

tickinfo <- map_dfr(tick_dates, function(d) {
  s <- combined_rpkmf %>% 
    filter(Date >= d) %>% 
    arrange(Date) %>% 
    slice(1)
  if(nrow(s) > 0) {
    tibble(sample_label = s$sample_label, label = format(d, "%b %Y"))
  } else {
    NULL
  }
})
tickinfo <- distinct(tickinfo, sample_label, .keep_all = TRUE)

# Define color palette
custom_colors <- c(
  "#FF7F00",  # Vivid orange
  "#00CED1",  # Dark turquoise
  "#BF6347",  # Tomato red
  "#6A5ACD",  # Slate blue
  "#B22222",  # Firebrick
  "#00BFFF",  # Deep sky blue
  "#C71585",  # Violet red
  "#3CB371",  # Medium sea green
  "#D2691E",  # Chocolate
  "#FF00FF",  # Magenta
  "#FF4500",  # Orange red
  "#00FF7F",  # Spring green
  "#7FFF00",  # Chartreuse
  "#FF6EB4",  # Pink
  "#484848",  # Graphite gray
  "#20B2AA",  # Light sea green
  "#FFDAB9",  # Peach puff
  "#800000",  # Maroon
  "#1E1E1E"   # Nearly black
)

if(length(top_strains) > length(custom_colors)) {
  stop("Too many strains for color palette")
}

strain_colors <- c(
  setNames(custom_colors[1:length(top_strains)], top_strains),
  "Other" = "#696969",
  "No Norwalk Detected" = "white"
)

# Create plots

# Panel 1: Log10 Norovirus read counts
panel1 <- ggplot(reads_per_sample, aes(x = sample_label, y = log10_reads)) +
  geom_col(fill = "#333333") +
  scale_y_continuous(
    breaks = 0:6,
    labels = as.character(0:6)
  ) +
  scale_x_discrete(
    breaks = tickinfo$sample_label,
    labels = NULL
  ) +
  labs(y = expression(Log[10]~"Norovirus Reads"), x = NULL, title = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.ticks.x.bottom = element_line(size = 0.8, color = "black"),
    axis.ticks.length.x = unit(0.18, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.text.y = element_text(size = 12),
    axis.ticks.y.left = element_line(size = 0.8, color = "black"),
    axis.line.y.left = element_line(color = "black", linewidth = 0.8),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

# Panel 2: Relative abundance of strains
panel2 <- ggplot(combined_rpkmf, 
                 aes(x = sample_label, y = relative_abundance, fill = strain_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = strain_colors, 
                    name = "Norwalk Virus Strain", 
                    drop = FALSE) +
  labs(
    x = "Sample ID and Date",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    labels = paste0(seq(0, 100, by = 20), "%"),
    expand = expansion(mult = c(0, 0.03)),
    limits = c(0, 100.5)
  ) +
  scale_x_discrete(
    breaks = tickinfo$sample_label,
    labels = tickinfo$label
  )

# Combine panels
final_plot <- plot_grid(
  panel1, panel2, 
  ncol = 1, 
  align = "v", 
  axis = "lr", 
  rel_heights = c(0.4, 1)
)

print(final_plot)


