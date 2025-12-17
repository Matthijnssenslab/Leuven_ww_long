# Rotavirus Genotype Analysis: Clinical vs Wastewater Surveillance
# Requires pre-processed data from prepare_rotavirus_data.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(grid)
library(ggpubr)
library(ggrepel)



# Load pre-processed datasets (in the local data/ folder)
wastewater_data <- readRDS("data/processed_wastewater_rotavirus.rds")
clinical_data <- readRDS("data/processed_clinical_rotavirus.rds")

# Define color schemes and order
g_level_order <- c("G1", "G2", "G8", "G9", 
                   "human-like G3", "equine-like G3", "other")

base_g_type_colors <- c(
  "G1" = "#FFD700",
  "G2" = "#228B22",
  "G8" = "#8FBC8F",
  "G9" = "#DC143C",
  "human-like G3" = "darkblue",
  "equine-like G3" = "#0976D2",
  "other" = "grey70"
)

base_p_type_colors <- c(
  "P[3]" = "#4169E1",
  "P[4]" = "#9370DB",
  "P[6]" = "#20B2AA",
  "P[8]" = "#FF1493",
  "other" = "grey70"
)

# Monthly aggregation
# Wastewater G-types
ww_g_month <- wastewater_data %>%
  filter(segment == "vp7") %>%
  group_by(month, geno_group) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  group_by(month) %>%
  mutate(month_total = sum(total_reads),
         percent = ifelse(month_total == 0, 0, total_reads / month_total * 100)) %>%
  ungroup() %>%
  mutate(geno_group = factor(geno_group, levels = g_level_order))

# Clinical G-types
cl_g_month <- clinical_data %>%
  filter(!is.na(g_type_harmonized)) %>%
  group_by(month, g_type_harmonized) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(month) %>%
  mutate(month_total = sum(count),
         percent = ifelse(month_total == 0, 0, count / month_total * 100)) %>%
  ungroup() %>%
  mutate(g_type_harmonized = factor(g_type_harmonized, levels = g_level_order))

# Wastewater P-types
ww_p_month <- wastewater_data %>%
  filter(segment == "vp4") %>%
  group_by(month, geno_group) %>%
  summarise(total_reads = sum(reads), .groups = "drop") %>%
  group_by(month) %>%
  mutate(month_total = sum(total_reads),
         percent = ifelse(month_total == 0, 0, total_reads / month_total * 100)) %>%
  ungroup()

# Clinical P-types
cl_p_month <- clinical_data %>%
  filter(!is.na(p_type_harmonized)) %>%
  group_by(month, p_type_harmonized) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(month) %>%
  mutate(month_total = sum(count),
         percent = ifelse(month_total == 0, 0, count / month_total * 100)) %>%
  ungroup()

# Complete missing month/genotype combinations
cl_g_month <- cl_g_month %>%
  complete(month, g_type_harmonized, 
           fill = list(count = 0, percent = 0, month_total = 0)) %>%
  group_by(month) %>%
  mutate(month_total = max(month_total)) %>%
  ungroup() %>%
  mutate(year = as.numeric(substr(as.character(month), 1, 4)))

cl_p_month <- cl_p_month %>%
  complete(month, p_type_harmonized, 
           fill = list(count = 0, percent = 0, month_total = 0)) %>%
  group_by(month) %>%
  mutate(month_total = max(month_total)) %>%
  ungroup() %>%
  mutate(year = as.numeric(substr(as.character(month), 1, 4)))

ww_g_month <- ww_g_month %>%
  complete(month, geno_group, 
           fill = list(total_reads = 0, percent = 0, month_total = 0)) %>%
  group_by(month) %>%
  mutate(month_total = max(month_total)) %>%
  ungroup() %>%
  mutate(year = as.numeric(substr(as.character(month), 1, 4)))

ww_p_month <- ww_p_month %>%
  complete(month, geno_group, 
           fill = list(total_reads = 0, percent = 0, month_total = 0)) %>%
  group_by(month) %>%
  mutate(month_total = max(month_total)) %>%
  ungroup() %>%
  mutate(year = as.numeric(substr(as.character(month), 1, 4)))

# Finalize color palettes
all_g_types_month <- levels(ww_g_month$geno_group)
all_p_types_month <- sort(unique(c(
  as.character(ww_p_month$geno_group), 
  as.character(cl_p_month$p_type_harmonized)
)))

final_g_type_colors_month <- setNames(
  sapply(all_g_types_month, function(gt) {
    if (gt %in% names(base_g_type_colors)) base_g_type_colors[gt] else "grey60"
  }),
  all_g_types_month
)

final_p_type_colors_month <- setNames(
  sapply(all_p_types_month, function(pt) {
    if (pt %in% names(base_p_type_colors)) base_p_type_colors[pt] else "grey60"
  }),
  all_p_types_month
)

# Create x-axis breaks (every 6 months)
all_months <- sort(unique(c(
  as.character(cl_g_month$month), 
  as.character(cl_p_month$month),
  as.character(ww_g_month$month), 
  as.character(ww_p_month$month)
)))

month_breaks <- all_months[seq(1, length(all_months), by = 6)]

# Align month factors
cl_g_month$month <- factor(cl_g_month$month, levels = all_months)
cl_p_month$month <- factor(cl_p_month$month, levels = all_months)
ww_g_month$month <- factor(ww_g_month$month, levels = all_months)
ww_p_month$month <- factor(ww_p_month$month, levels = all_months)

# Publication themes
theme_pub <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(2, "pt"),
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black", size = 0.3),
      axis.line.x = element_line(color = "black", size = 0.3),
      plot.margin = unit(c(2, 2, 2, 2), "pt")
    )
}

theme_pub_box <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.ticks.length = unit(2, "pt"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3),
      plot.margin = unit(c(2, 2, 2, 2), "pt")
    )
}

# Calculate axis limits
max_cl_cases <- max(
  max(cl_g_month$month_total, na.rm = TRUE), 
  max(cl_p_month$month_total, na.rm = TRUE)
) * 1.1

max_ww_reads_log <- max(
  log10(max(ww_g_month$month_total, na.rm = TRUE) + 1),
  log10(max(ww_p_month$month_total, na.rm = TRUE) + 1)
) * 1.1

# Create time series plots
# G-type plots
p_g_cl_count_pub <- ggplot(
  cl_g_month %>% distinct(month, month_total),
  aes(x = month, y = month_total)
) +
  geom_col(fill = "black", width = 0.8) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  coord_cartesian(ylim = c(0, max_cl_cases), expand = FALSE) +
  theme_pub()

p_g_cl_rel_pub <- ggplot(
  cl_g_month, 
  aes(x = month, y = percent, fill = g_type_harmonized)
) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = final_g_type_colors_month, na.translate = FALSE) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), oob = scales::squish) +
  guides(fill = "none") +
  theme_pub_box()

p_g_ww_count_pub <- ggplot(
  ww_g_month %>% 
    distinct(month, month_total) %>% 
    mutate(reads_log = log10(month_total + 1)),
  aes(x = month, y = reads_log)
) +
  geom_col(fill = "black", width = 0.8) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  coord_cartesian(ylim = c(0, max_ww_reads_log), expand = FALSE) +
  theme_pub()

p_g_ww_rel_pub <- ggplot(
  ww_g_month, 
  aes(x = month, y = percent, fill = geno_group)
) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = final_g_type_colors_month, na.translate = FALSE) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), oob = scales::squish) +
  guides(fill = "none") +
  theme_pub_box()

# P-type plots
p_p_cl_count_pub <- ggplot(
  cl_p_month %>% distinct(month, month_total),
  aes(x = month, y = month_total)
) +
  geom_col(fill = "black", width = 0.8) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  coord_cartesian(ylim = c(0, max_cl_cases), expand = FALSE) +
  theme_pub()

p_p_cl_rel_pub <- ggplot(
  cl_p_month, 
  aes(x = month, y = percent, fill = p_type_harmonized)
) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = final_p_type_colors_month, na.translate = FALSE) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), oob = scales::squish) +
  guides(fill = "none") +
  theme_pub_box()

p_p_ww_count_pub <- ggplot(
  ww_p_month %>% 
    distinct(month, month_total) %>% 
    mutate(reads_log = log10(month_total + 1)),
  aes(x = month, y = reads_log)
) +
  geom_col(fill = "black", width = 0.8) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  coord_cartesian(ylim = c(0, max_ww_reads_log), expand = FALSE) +
  theme_pub()

p_p_ww_rel_pub <- ggplot(
  ww_p_month, 
  aes(x = month, y = percent, fill = geno_group)
) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = final_p_type_colors_month, na.translate = FALSE) +
  scale_x_discrete(expand = c(0, 0), breaks = month_breaks) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), oob = scales::squish) +
  guides(fill = "none") +
  theme_pub_box()

# Combine time series figure

fig_publication <- grid.arrange(
  arrangeGrob(
    p_g_cl_count_pub,
    p_g_cl_rel_pub,
    p_g_ww_count_pub,
    p_g_ww_rel_pub,
    ncol = 1, heights = c(0.5, 1.5, 0.5, 1.5)
  ),
  arrangeGrob(
    p_p_cl_count_pub,
    p_p_cl_rel_pub,
    p_p_ww_count_pub,
    p_p_ww_rel_pub,
    ncol = 1, heights = c(0.5, 1.5, 0.5, 1.5)
  ),
  ncol = 2, widths = c(1, 1)
)

# Prepare lollipop plot data
# Define season months (>50 clinical cases)
season_months <- cl_g_month %>%
  filter(month_total > 50) %>%
  distinct(year, month)

# Wastewater G-types: per-sample relative abundance
ww_g_per_sample <- wastewater_data %>%
  filter(segment == "vp7") %>%
  group_by(sample_ID, month, year, geno_group) %>%
  summarise(sample_reads = sum(reads), .groups = "drop") %>%
  group_by(sample_ID, month, year) %>%
  mutate(sample_total = sum(sample_reads),
         sample_percent = ifelse(sample_total == 0, 0, 
                                 sample_reads / sample_total * 100)) %>%
  ungroup() %>%
  group_by(month, year) %>%
  mutate(total_samples = n_distinct(sample_ID)) %>%
  group_by(month, year, geno_group, total_samples) %>%
  summarise(sum_percent = sum(sample_percent), .groups = "drop") %>%
  mutate(avg_percent = sum_percent / total_samples) %>%
  select(month, year, geno_group, avg_percent) %>%
  filter(geno_group != "other")

# Wastewater P-types: per-sample relative abundance
ww_p_per_sample <- wastewater_data %>%
  filter(segment == "vp4") %>%
  group_by(sample_ID, month, year, geno_group) %>%
  summarise(sample_reads = sum(reads), .groups = "drop") %>%
  group_by(sample_ID, month, year) %>%
  mutate(sample_total = sum(sample_reads),
         sample_percent = ifelse(sample_total == 0, 0, 
                                 sample_reads / sample_total * 100)) %>%
  ungroup() %>%
  mutate(geno_group = ifelse(geno_group %in% c("P[6]", "P[3]"), 
                             "other", geno_group)) %>%
  group_by(month, year) %>%
  mutate(total_samples = n_distinct(sample_ID)) %>%
  group_by(month, year, geno_group, total_samples) %>%
  summarise(sum_percent = sum(sample_percent), .groups = "drop") %>%
  mutate(avg_percent = sum_percent / total_samples) %>%
  select(month, year, geno_group, avg_percent) %>%
  filter(geno_group %in% c(cl_p_month$p_type_harmonized, "other") | 
           geno_group == "other") %>%
  filter(geno_group != "other")

# Clinical seasonal data
cl_g_season <- cl_g_month %>%
  inner_join(season_months, by = c("year", "month")) %>%
  filter(g_type_harmonized != "other") %>%
  select(month, year, genotype = g_type_harmonized, avg_percent = percent)

cl_p_season <- cl_p_month %>%
  inner_join(season_months, by = c("year", "month")) %>%
  filter(p_type_harmonized != "other") %>%
  select(month, year, genotype = p_type_harmonized, avg_percent = percent)

# Wastewater seasonal data
ww_g_season <- ww_g_per_sample %>%
  semi_join(season_months, by = c("year", "month")) %>%
  select(month, year, genotype = geno_group, avg_percent)

ww_p_season <- ww_p_per_sample %>%
  semi_join(season_months, by = c("year", "month")) %>%
  select(month, year, genotype = geno_group, avg_percent)

# Merge and order lollipop data
lollipop_data <- bind_rows(
  cl_g_season %>% 
    group_by(year, genotype) %>%
    summarise(avg_percent = mean(avg_percent), .groups = "drop") %>%
    mutate(type = "G", source = "Clinical", x_offset = -0.16),
  
  ww_g_season %>% 
    group_by(year, genotype) %>%
    summarise(avg_percent = mean(avg_percent), .groups = "drop") %>%
    mutate(type = "G", source = "Wastewater", x_offset = -0.05),
  
  cl_p_season %>% 
    group_by(year, genotype) %>%
    summarise(avg_percent = mean(avg_percent), .groups = "drop") %>%
    mutate(type = "P", source = "Clinical", x_offset = 0.05),
  
  ww_p_season %>% 
    group_by(year, genotype) %>%
    summarise(avg_percent = mean(avg_percent), .groups = "drop") %>%
    mutate(type = "P", source = "Wastewater", x_offset = 0.16)
) %>% 
  filter(avg_percent > 0) %>%
  mutate(
    year_label = paste("Year", year),
    genotype = as.character(genotype),
    source = factor(source, levels = c("Clinical", "Wastewater"))
  )

# Order by type (G then P), then by wastewater highest %
ww_max <- lollipop_data %>%
  filter(source == "Wastewater") %>%
  group_by(type, genotype) %>%
  summarise(max_ww = max(avg_percent), .groups = "drop") %>%
  arrange(type, desc(max_ww))

order_df <- ww_max %>% mutate(geno_id = paste(type, genotype))

lollipop_data <- lollipop_data %>%
  mutate(geno_id = paste(type, genotype)) %>%
  mutate(geno_id = factor(geno_id, levels = order_df$geno_id))

# Create lollipop plot
fig_lollipop_pub <- ggplot(
  lollipop_data,
  aes(x = as.numeric(geno_id) + x_offset, y = avg_percent)
) +
  geom_segment(
    aes(xend = as.numeric(geno_id) + x_offset, yend = 0, color = source),
    size = 1, alpha = 0.5
  ) +
  geom_point(
    aes(color = source, fill = source, shape = source),
    size = 4, stroke = 1, alpha = 0.9
  ) +
  geom_text_repel(
    aes(label = round(avg_percent, 1), color = source),
    size = 2.3, fontface = "bold", show.legend = FALSE,
    direction = "both", seed = 42, max.overlaps = 100,
    box.padding = 0.35, segment.color = "grey60"
  ) +
  scale_color_manual(
    values = c("Clinical" = "#D4AF37", "Wastewater" = "#6B2C91"), 
    name = "Source"
  ) +
  scale_fill_manual(
    values = c("Clinical" = "#D4AF37", "Wastewater" = "#6B2C91"), 
    name = "Source"
  ) +
  scale_shape_manual(
    values = c("Clinical" = 21, "Wastewater" = 22), 
    name = "Source"
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(lollipop_data$geno_id)),
    labels = gsub("^[A-Z] ", "", levels(lollipop_data$geno_id)),
    expand = expansion(mult = c(0.09, 0.16))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.13)), 
    limits = c(0, 100)
  ) +
  facet_wrap(~year_label, ncol = 2) +
  coord_flip() +
  labs(
    title = "Seasonal Peak Genotypes (G-types & P-types): Clinical vs Wastewater",
    x = "", 
    y = "Average RA (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                              margin = margin(b = 10)),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 5)),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.x = unit(3, "pt"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length.y = unit(3, "pt"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    strip.background = element_rect(fill = "gray", color = "black", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Display plots
print(fig_publication)
print(fig_lollipop_pub)


