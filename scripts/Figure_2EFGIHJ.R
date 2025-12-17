## Mastadenovirus: strain and species relative abundance vs clinical cases

# Setup
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  library(scales)
  library(zoo)
  library(ggpubr)
  library(patchwork)
})

# Wastewater data (pre-filtered dataset)
filtered_data <- readRDS("data/filtered_data_wwlong.rds")

all_samples <- filtered_data %>% 
  distinct(sample_ID, Date)

# ------------------------------------------------------------------
# 2. Prepare Mastadenovirus data
#BLASTN shows this unclassified strain consensus sequence as HAdV-41 therefore we interpret it as HAdV-41
mastadenovirus_data <- filtered_data %>% 
  filter(genus == "Mastadenovirus") %>%
  mutate(strain = if_else(
    strain == "unclassified Human mastadenovirus F strain",
    "Human adenovirus 41",
    strain
  ))

#We remove adenovirus strains
species_to_remove <- c("Murine mastadenovirus B", 
                       "Porcine mastadenovirus B", 
                       "Porcine mastadenovirus C")

mastadenovirus_data <- mastadenovirus_data %>%
  filter(!species %in% species_to_remove)

# ------------------------------------------------------------------
# 3. Aggregate RPKMF within Mastadenovirus
# ------------------------------------------------------------------
strain_rpkmf <- mastadenovirus_data %>%
  group_by(sample_ID, Date, strain) %>%
  summarise(total_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

no_detection_samples <- all_samples %>%
  anti_join(strain_rpkmf, by = c("sample_ID", "Date")) %>%
  mutate(strain = "No Mastadenovirus Detected", total_RPKMF = 0)

combined_rpkmf <- bind_rows(strain_rpkmf, no_detection_samples)

# ------------------------------------------------------------------
# 4. Assign top strains and "Other"
# ------------------------------------------------------------------
top_n <- 10

top_strains <- combined_rpkmf %>%
  group_by(strain) %>%
  summarise(total = sum(total_RPKMF, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  filter(strain != "No Mastadenovirus Detected") %>%
  slice(1:top_n) %>%
  pull(strain)

custom_strain_order <- c("Human adenovirus 41")
other_strains <- setdiff(top_strains, custom_strain_order)
final_strain_order <- c(custom_strain_order, other_strains, "Other")

combined_rpkmf <- combined_rpkmf %>%
  mutate(strain_grouped = case_when(
    strain == "No Mastadenovirus Detected" ~ "No Mastadenovirus Detected",
    strain %in% custom_strain_order ~ strain,
    strain %in% other_strains ~ strain,
    TRUE ~ "Other"
  )) %>%
  group_by(sample_ID, Date, strain_grouped) %>%
  summarise(RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

combined_rpkmf$strain_grouped <- factor(
  combined_rpkmf$strain_grouped,
  levels = final_strain_order
)

# ------------------------------------------------------------------
# 5. Date-based aggregation
# ------------------------------------------------------------------
date_aggregated <- combined_rpkmf %>%
  group_by(Date, strain_grouped) %>%
  summarise(
    avg_RPKMF = mean(RPKMF, na.rm = TRUE),
    total_samples = n_distinct(sample_ID),
    .groups = "drop"
  ) %>%
  group_by(Date) %>%
  mutate(
    total_date_RPKMF = sum(avg_RPKMF, na.rm = TRUE),
    relative_abundance = ifelse(total_date_RPKMF > 0, 
                                (avg_RPKMF / total_date_RPKMF) * 100, 0)
  ) %>%
  ungroup()

date_aggregated$strain_grouped <- factor(
  date_aggregated$strain_grouped,
  levels = final_strain_order
)

# ------------------------------------------------------------------
# 6. Load clinical data and aggregate
#    Expecting: data/public_cases.csv (epistat export, aggregated counts)
# ------------------------------------------------------------------
epistat <- read_csv("data/public_cases.csv") %>%
  filter(Subject == "V_ADV", 
         DateMonday >= as.Date("2022-04-24"), 
         DateMonday <= as.Date("2024-12-13"))

epistat_weekly <- epistat %>%
  mutate(Week = isoweek(DateMonday), Year = year(DateMonday)) %>%
  group_by(Year, Week) %>%
  summarise(Date = min(DateMonday), cases = n(), .groups = "drop")

wastewater_dates <- combined_rpkmf %>% 
  distinct(Date) %>% 
  arrange(Date)

cases_by_date <- epistat_weekly %>%
  mutate(Date = as.Date(sapply(Date, function(d) {
    wastewater_dates$Date[which.min(abs(wastewater_dates$Date - d))]
  }))) %>%
  group_by(Date) %>%
  summarise(cases = mean(cases), .groups = "drop")

# ------------------------------------------------------------------
# 7. Date index and x-axis labels
# ------------------------------------------------------------------
date_order <- date_aggregated %>%
  distinct(Date) %>%
  arrange(Date) %>%
  mutate(date_index = row_number())

date_aggregated <- date_aggregated %>%
  left_join(date_order, by = "Date")

cases_by_date <- cases_by_date %>%
  left_join(date_order, by = "Date")

start_date <- floor_date(min(date_order$Date), "month")
end_date <- ceiling_date(max(date_order$Date), "month") - days(1)

bi_month_seq <- seq.Date(start_date, end_date, by = "2 months")
bi_month_breaks <- sapply(bi_month_seq, function(d) {
  s <- date_order %>% filter(Date >= d) %>% slice(1)
  if (nrow(s) > 0) s$date_index else NA
})
bi_month_breaks <- bi_month_breaks[!is.na(bi_month_breaks)]

bi_month_labels <- sapply(bi_month_breaks, function(idx) {
  date_val <- date_order %>% filter(date_index == idx) %>% pull(Date)
  format(date_val, "%b-%Y")
})

# ------------------------------------------------------------------
# 8. Plot colors
# ------------------------------------------------------------------
palette <- c("#a8ceee", "#E41A1C", "#FF7F00", "#984EA3", "#4DAF4A",
             "#FFFF33", "#A65628", "#F781BF", "#999999", "#fb8072", 
             "black", "white")
names(palette) <- final_strain_order
strain_colors <- palette

# ------------------------------------------------------------------
# 9. Dual-axis scaling
# ------------------------------------------------------------------
scale_factor <- ifelse(max(cases_by_date$cases, na.rm = TRUE) > 0,
                       max(cases_by_date$cases, na.rm = TRUE) / 100, 1)

# ------------------------------------------------------------------
# 10. Panel A: Strain-level relative abundance
# ------------------------------------------------------------------
strains_relative <- ggplot() +
  geom_col(
    data = date_aggregated,
    aes(x = date_index, y = relative_abundance, fill = strain_grouped),
    width = 1
  ) +
  scale_fill_manual(
    values = strain_colors,
    limits = final_strain_order,
    drop = FALSE,
    name = "Mastadenovirus strain"
  ) +
  geom_line(
    data = cases_by_date,
    aes(x = date_index, y = cases / scale_factor),
    color = "black", size = 1
  ) +
  geom_point(
    data = cases_by_date,
    aes(x = date_index, y = cases / scale_factor),
    color = "black", size = 2
  ) +
  scale_x_continuous(
    name = "Date",
    breaks = bi_month_breaks,
    labels = bi_month_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Relative Abundance (%)",
    labels = percent_format(scale = 1),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * scale_factor, name = "Case Numbers")
  ) +
  labs(
    title = "Mastadenovirus: Relative Abundance & Clinical Cases Over Time",
    subtitle = paste("Bars:", paste(top_strains, collapse = ", "), 
                     ", Other | Black line: clinical cases")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )

strains_relative

# ------------------------------------------------------------------
# 11. Species-level aggregation
# ------------------------------------------------------------------
species_rpkmf <- mastadenovirus_data %>%
  group_by(sample_ID, Date, species) %>%
  summarise(total_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

no_detection_samples_species <- all_samples %>%
  anti_join(species_rpkmf, by = c("sample_ID", "Date")) %>%
  mutate(species = "No Mastadenovirus Detected", total_RPKMF = 0)

combined_rpkmf_species <- bind_rows(species_rpkmf, no_detection_samples_species)

species_totals <- combined_rpkmf_species %>%
  group_by(species) %>%
  summarise(total = sum(total_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  filter(species != "No Mastadenovirus Detected")

custom_species_order <- "Human mastadenovirus F"
other_species <- setdiff(species_totals$species, custom_species_order)
final_species_order <- c(custom_species_order, other_species, "Other")

combined_rpkmf_species <- combined_rpkmf_species %>%
  mutate(species_grouped = case_when(
    species == "No Mastadenovirus Detected" ~ "No Mastadenovirus Detected",
    species %in% custom_species_order ~ species,
    species %in% other_species ~ species,
    TRUE ~ "Other"
  )) %>%
  group_by(sample_ID, Date, species_grouped) %>%
  summarise(RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

combined_rpkmf_species$species_grouped <- factor(
  combined_rpkmf_species$species_grouped,
  levels = final_species_order
)

date_aggregated_species <- combined_rpkmf_species %>%
  group_by(Date, species_grouped) %>%
  summarise(
    avg_RPKMF = mean(RPKMF, na.rm = TRUE),
    total_samples = n_distinct(sample_ID),
    .groups = "drop"
  ) %>%
  group_by(Date) %>%
  mutate(
    total_date_RPKMF = sum(avg_RPKMF, na.rm = TRUE),
    relative_abundance = ifelse(total_date_RPKMF > 0, 
                                (avg_RPKMF / total_date_RPKMF) * 100, 0)
  ) %>%
  ungroup()

date_aggregated_species$species_grouped <- factor(
  date_aggregated_species$species_grouped,
  levels = final_species_order
)

date_aggregated_species <- date_aggregated_species %>%
  left_join(date_order, by = "Date")

pastel_palette <- c(
  "#9ecae6", "#f9caf4", "#4DAF4A", "gray", "#264653", 
  "#c2fae0", "#c4d96f", "#b4b8ab", "#e9c46a", "#264653", 
  "#22223b", "#264653"
)

if(length(final_species_order) > length(pastel_palette)) {
  pastel_palette <- rep(pastel_palette, length.out = length(final_species_order))
}
names(pastel_palette) <- final_species_order
species_colors <- pastel_palette[final_species_order]

# ------------------------------------------------------------------
# 12. Panel B: Species-level relative abundance
# ------------------------------------------------------------------
species_relative <- ggplot() +
  geom_col(
    data = date_aggregated_species,
    aes(x = date_index, y = relative_abundance, fill = species_grouped),
    width = 1
  ) +
  scale_fill_manual(
    values = species_colors,
    limits = final_species_order,
    drop = FALSE,
    name = "Mastadenovirus species"
  ) +
  geom_line(
    data = cases_by_date,
    aes(x = date_index, y = cases / scale_factor),
    color = "#14213d", size = 1
  ) +
  geom_point(
    data = cases_by_date,
    aes(x = date_index, y = cases / scale_factor),
    color = "#14213d", size = 2
  ) +
  scale_x_continuous(
    name = "Date",
    breaks = bi_month_breaks,
    labels = bi_month_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Relative Abundance (%)",
    labels = percent_format(scale = 1),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * scale_factor, name = "Case Numbers")
  ) +
  labs(
    title = "Mastadenovirus: Species Relative Abundance & Clinical Cases Over Time",
    subtitle = paste("Bars:", paste(final_species_order, collapse = ", "), 
                     "| Line: clinical cases")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

species_relative

# ------------------------------------------------------------------
# 13. Correlation calculations
# ------------------------------------------------------------------
priority_strains <- custom_strain_order

sample_sums <- combined_rpkmf %>%
  group_by(sample_ID, Date) %>%
  summarise(
    total_adeno_RPKMF = sum(RPKMF[!strain_grouped %in% c("No Mastadenovirus Detected")], 
                            na.rm = TRUE),
    resp_adeno_RPKMF = sum(RPKMF[!(strain_grouped %in% c(priority_strains, 
                                                         "No Mastadenovirus Detected", 
                                                         "Other"))], 
                           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    resp_adeno_pct = ifelse(total_adeno_RPKMF > 0, 
                            (resp_adeno_RPKMF / total_adeno_RPKMF) * 100, 0)
  )

resp_fraction_by_date <- sample_sums %>%
  group_by(Date) %>%
  summarise(resp_adeno_pct = mean(resp_adeno_pct, na.rm = TRUE), .groups = "drop")

corr_df <- resp_fraction_by_date %>%
  left_join(cases_by_date, by = "Date") %>%
  filter(!is.na(resp_adeno_pct), !is.na(cases))

# ------------------------------------------------------------------
# 14. Correlation scatter plots by window
# ------------------------------------------------------------------
roll_mean_na <- function(x, k = 4, partial = TRUE) {
  rollapply(x, width = k, 
            FUN = function(z) mean(z, na.rm = TRUE), 
            fill = NA, align = "right", partial = partial)
}

make_window_plot <- function(df, title) {
  ggplot(df, aes(x = resp_2wk, y = cases_2wk)) +
    geom_point(color = "darkblue", size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, 
                color = "orange", fill = "orange", alpha = 0.1) +
    stat_cor(method = "spearman", size = 5, label.x = 5) +
    labs(
      x = "Respiratory Fraction (2-wk avg)", 
      y = "Clinical Cases (2-wk avg)", 
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.25, "cm")
    )
}

dfA <- corr_df %>%
  filter(Date >= as.Date("2022-04-25"), Date <= as.Date("2022-08-01")) %>%
  arrange(Date) %>%
  mutate(
    resp_2wk = roll_mean_na(resp_adeno_pct, k = 4, partial = TRUE),
    cases_2wk = roll_mean_na(cases, k = 4, partial = TRUE)
  ) %>%
  filter(!is.na(resp_2wk), !is.na(cases_2wk))

pA <- make_window_plot(dfA, "Window A: Apr–Aug 2022")

dfB <- corr_df %>%
  filter(Date >= as.Date("2022-09-15"), Date <= as.Date("2023-08-03")) %>%
  arrange(Date) %>%
  mutate(
    resp_2wk = roll_mean_na(resp_adeno_pct, k = 4, partial = TRUE),
    cases_2wk = roll_mean_na(cases, k = 4, partial = TRUE)
  ) %>%
  filter(!is.na(resp_2wk), !is.na(cases_2wk))

pB <- make_window_plot(dfB, "Window B: Sep 2022–Aug 2023")

dfC <- corr_df %>%
  filter(Date >= as.Date("2023-10-01"), Date <= as.Date("2024-05-01")) %>%
  arrange(Date) %>%
  mutate(
    resp_2wk = roll_mean_na(resp_adeno_pct, k = 4, partial = TRUE),
    cases_2wk = roll_mean_na(cases, k = 4, partial = TRUE)
  ) %>%
  filter(!is.na(resp_2wk), !is.na(cases_2wk))

pC <- make_window_plot(dfC, "Window C: Oct 2023–May 2024")

dfD <- corr_df %>%
  filter(Date >= as.Date("2024-06-01"), Date <= as.Date("2024-12-13")) %>%
  arrange(Date) %>%
  mutate(
    resp_2wk = roll_mean_na(resp_adeno_pct, k = 4, partial = TRUE),
    cases_2wk = roll_mean_na(cases, k = 4, partial = TRUE)
  ) %>%
  filter(!is.na(resp_2wk), !is.na(cases_2wk))

pD <- make_window_plot(dfD, "Window D: Jun–Dec 2024")

combined_plot <- pA + pB + pC + pD + plot_layout(ncol = 4, nrow = 1)
combined_plot

# ------------------------------------------------------------------
# 15. Window annotations for time series (using nearest date matching)
# ------------------------------------------------------------------
find_closest_index <- function(target_date, date_df) {
  date_df %>%
    mutate(diff = abs(Date - target_date)) %>%
    arrange(diff) %>%
    slice(1) %>%
    pull(date_index)
}

window_dates <- tibble(
  win_id = c("A", "B", "C", "D"),
  date_start = as.Date(c("2022-04-25", "2022-09-15", "2023-10-01", "2024-06-01")),
  date_end = as.Date(c("2022-08-01", "2023-08-03", "2024-05-01", "2024-12-13"))
)

window_bounds <- window_dates %>%
  rowwise() %>%
  mutate(
    idx_start = find_closest_index(date_start, date_order),
    idx_end = find_closest_index(date_end, date_order)
  ) %>%
  ungroup()

window_colors <- c("A" = "#E41A1C", "B" = "#377EB8", 
                   "C" = "#4DAF4A", "D" = "#984EA3")

half_index <- 0.5

window_vlines <- c(
  lapply(window_bounds$idx_start, function(idx)
    geom_vline(xintercept = idx, linetype = "dashed", 
               color = "black", linewidth = 1)),
  lapply(window_bounds$idx_end, function(idx)
    geom_vline(xintercept = idx, linetype = "dashed", 
               color = "black", linewidth = 1))
)

species_relative_annotated <- species_relative +
  window_vlines +
  labs(title = NULL, subtitle = NULL) +
  theme(plot.margin = margin(t = -5, r = 10, b = 25, l = 10))

strains_relative_no_x_window <- strains_relative +
  window_vlines +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 0, r = 10, b = -5, l = 10)
  )

strip_data <- data.frame(
  xstart = window_bounds$idx_start - half_index,
  xend   = window_bounds$idx_end + half_index,
  win_id = window_bounds$win_id
)

window_strip <- ggplot() +
  geom_rect(
    data = strip_data,
    aes(xmin = xstart, xmax = xend, ymin = 0, ymax = 1, fill = win_id),
    color = NA, alpha = 1
  ) +
  scale_fill_manual(values = setNames(window_colors, window_bounds$win_id), 
                    name = "Window") +
  scale_x_continuous(
    breaks = bi_month_breaks, 
    labels = bi_month_labels, 
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = NULL, 
    labels = NULL, 
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

figure_timeseries_with_windows <-
  strains_relative_no_x_window /
  species_relative_annotated /
  window_strip +
  plot_layout(heights = c(2, 2, 0.08), guides = "collect") &
  theme(legend.position = "right")

figure_timeseries_with_windows

figure_correlation_panel <- combined_plot
figure_correlation_panel
