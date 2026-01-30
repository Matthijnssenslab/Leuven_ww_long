# Figure 1E: Relative Abundance of Viral Families Over Time (Belgium)

# ------------------------------
# 0. Setup
# ------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# 1. Data loading
filtered_data <- readRDS("data/filtered_data_wwlong.rds")

# ------------------------------
family_reads <- filtered_data %>%
  group_by(sample_ID, Date, family) %>%
  summarise(
    reads_mapped = sum(reads_aligned, na.rm = TRUE),
    .groups = "drop"
  )

total_reads <- family_reads %>%
  group_by(sample_ID, Date) %>%
  summarise(
    total_mapped_reads = sum(reads_mapped, na.rm = TRUE),
    .groups = "drop"
  )

relative_abundance <- family_reads %>%
  left_join(total_reads, by = c("sample_ID", "Date")) %>%
  mutate(
    relative_abundance = 100 * reads_mapped / pmax(total_mapped_reads, 1)
  )

# Establish a stable sample order (by Date, then sample_ID)
order_df <- relative_abundance %>%
  distinct(sample_ID, Date) %>%
  arrange(Date, sample_ID) %>%
  mutate(
    sample_order     = row_number(),
    sample_order_fct = factor(sample_order, levels = sample_order)
  )

# 3. Top families and grouping
top_families <- relative_abundance %>%
  group_by(family) %>%
  summarise(
    total_reads = sum(reads_mapped, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = 10) %>%
  pull(family)

plot_data <- relative_abundance %>%
  mutate(
    family_grouped = if_else(family %in% top_families, family, "Other")
  ) %>%
  group_by(sample_ID, Date, family_grouped) %>%
  summarise(
    relative_abundance = sum(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(order_df, by = c("sample_ID", "Date")) %>%
  arrange(sample_order)

# Month tick positions (first sample of each month)
month_labels_df <- order_df %>%
  mutate(ym = format(Date, "%Y-%m")) %>%
  group_by(ym) %>%
  summarise(
    idx        = first(sample_order),
    label_date = first(Date),
    .groups    = "drop"
  )

# 4. Palette (distinct, readable)
color_palette <- c(
  "Adenoviridae"     = "#E41A1C",
  "Astroviridae"     = "#377EB8",
  "Caliciviridae"    = "#4DAF4A",
  "Circoviridae"     = "#984EA3",
  "Dicistroviridae"  = "#FF7F00",
  "Herpesviridae"    = "#FFFF33",
  "Papillomaviridae" = "#A65628",
  "Parvoviridae"     = "#F781BF",
  "Picornaviridae"   = "#999999",
  "Polyomaviridae"   = "#66C2A5",
  "Retroviridae"     = "#FC8D62",
  "Sedoreoviridae"   = "#8DA0CB",
  "Other"            = "black"
)

missing_cols <- setdiff(unique(plot_data$family_grouped), names(color_palette))
if (length(missing_cols) > 0) {
  extra_cols <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(missing_cols)),
    missing_cols
  )
  color_palette <- c(color_palette, extra_cols)
}

# 5. Bottom panel: stacked relative abundance
bottom_plot <- ggplot(
  plot_data,
  aes(x = factor(sample_order, levels = order_df$sample_order),
      y = relative_abundance,
      fill = family_grouped)
) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_x_discrete(
    breaks = as.character(month_labels_df$idx),
    labels = format(month_labels_df$label_date, "%b-%Y")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = color_palette, drop = FALSE) +
  labs(
    x    = "Sample date",
    y    = "Relative abundance (%)",
    fill = "Viral family"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, size = 0.8),
    axis.ticks         = element_line(color = "black"),
    axis.ticks.length  = unit(0.05, "cm"),
    legend.position    = "right",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9)
  )

# 6. Top panel: total reads per sample
total_reads_top <- filtered_data %>%
  group_by(sample_ID, Date) %>%
  summarise(
    total_reads_aligned = sum(reads_aligned, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(order_df, by = c("sample_ID", "Date")) %>%
  arrange(sample_order)

top_plot <- ggplot(
  total_reads_top,
  aes(x = factor(sample_order, levels = order_df$sample_order),
      y = total_reads_aligned)
) +
  geom_col(fill = "black", color = "black", width = 0.9, linewidth = 0.2) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = function(x) scales::number(x / 1e6, accuracy = 0.1)
  ) +
  labs(
    x = NULL,
    y = "Reads mapped (millions)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid          = element_blank(),
    panel.border        = element_blank(),
    axis.line.y.left    = element_line(color = "black"),
    axis.line.x.bottom  = element_line(color = "black"),
    axis.title.x        = element_blank(),
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    axis.ticks.y        = element_line(color = "black"),
    axis.ticks.length.y = unit(3, "pt"),
    axis.text.y         = element_text(color = "black"),
    plot.margin         = margin(b = 2, t = 2, r = 5, l = 8)
  )

# 7. Combine and save
fig <- top_plot / bottom_plot + plot_layout(heights = c(1, 4))

print(fig)
