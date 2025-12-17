# Dynamic Time Warping Analysis: Rotavirus and Adenovirus Lag Analysis
# Compares wastewater viral signals with clinical case data using DTW

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(lubridate)
library(fuzzyjoin)
library(dtw)
library(zoo)
library(patchwork)


# Load pre-filtered data
filtered_data <- readRDS("data/filtered_data_wwlong.rds")


# PART 1: ROTAVIRUS LAG ANALYSIS


# Load clinical data for Rotavirus
epistat <- read_csv("data/public_cases.csv") %>%
  filter(Subject == "V_RTV",
         DateMonday >= as.Date("2022-04-24") & 
           DateMonday <= as.Date("2024-12-17"))

# Aggregate clinical data by week
clinical_data_rota <- epistat %>%
  group_by(DateMonday) %>%
  summarise(clinical_cases = n(), .groups = "drop")

# Aggregate wastewater rotavirus RPKMF by date
wastewater_rotavirus <- filtered_data %>%
  filter(genus == "Rotavirus") %>%
  group_by(Date) %>%
  summarise(RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

# Fuzzy join clinical and wastewater data
joined_data_rota <- difference_left_join(
  clinical_data_rota, 
  wastewater_rotavirus,
  by = c("DateMonday" = "Date"),
  max_dist = 5,
  distance_col = "date_diff"
) %>%
  group_by(DateMonday) %>%
  filter(abs(date_diff) == min(abs(date_diff))) %>%
  ungroup() %>%
  arrange(DateMonday) %>%
  drop_na()

# Standardize and smooth data
joined_data_rota <- joined_data_rota %>%
  mutate(
    log_clinical = scale(log2(clinical_cases + 1)),
    log_rpkmf = scale(log10(RPKMF + 1))
  ) %>%
  mutate(
    smooth_log_clinical = rollmean(log_clinical, k = 2, fill = NA, align = "right"),
    smooth_log_rpkmf = rollmean(log_rpkmf, k = 2, fill = NA, align = "right")
  ) %>%
  drop_na()

# Perform DTW alignment
ts_wastewater_rota <- joined_data_rota$smooth_log_rpkmf
ts_clinical_rota <- joined_data_rota$smooth_log_clinical

dtw_alignment_rota <- dtw(
  ts_wastewater_rota, 
  ts_clinical_rota,
  keep = TRUE,
  step.pattern = symmetric2,
  window.type = "sakoechiba",
  window.size = 4
)

# Extract lags
lags_dtw_rota <- dtw_alignment_rota$index2 - dtw_alignment_rota$index1
lag_df_rota <- data.frame(lag = lags_dtw_rota)

# Create complete frequency table for lags
lag_counts_rota <- lag_df_rota %>%
  count(lag) %>%
  complete(lag = -4:4, fill = list(n = 0))

# Plot rotavirus lags
rotavirus_lag_plot <- ggplot(lag_counts_rota, aes(x = lag, y = n)) +
  geom_col(fill = "#345", color = "black", alpha = 0.4, width = 0.9) +
  scale_x_continuous(breaks = -4:4, limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(
    title = "Rotavirus",
    subtitle = "Negative lag: Clinical earlier; Positive lag: Wastewater leads",
    x = "Lag (clinical - wastewater index)",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  coord_cartesian(clip = "off")

# PART 2: ADENOVIRUS LAG ANALYSIS (3 WINDOWS)


cat("Analyzing Adenovirus lag across time windows...\n")

# Focus on Mastadenovirus
adenovirus_data <- filtered_data %>%
  filter(genus == "Mastadenovirus") %>%
  mutate(strain = if_else(
    strain == "unclassified Human mastadenovirus F strain",
    "Human adenovirus 41",
    strain
  ))

# Aggregate RPKMF and calculate relative abundances
strain_rpkmf <- adenovirus_data %>%
  group_by(sample_ID, Date, strain) %>%
  summarise(total_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

all_samples <- filtered_data %>% 
  distinct(sample_ID, Date)

no_detection_samples <- all_samples %>%
  anti_join(strain_rpkmf, by = c("sample_ID", "Date")) %>%
  mutate(strain = "No Mastadenovirus Detected", total_RPKMF = 0)

combined_rpkmf <- bind_rows(strain_rpkmf, no_detection_samples)

num_top_strains <- 12
top_strains <- combined_rpkmf %>%
  filter(strain != "No Mastadenovirus Detected") %>%
  group_by(strain) %>%
  summarise(total_RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_RPKMF)) %>%
  slice_head(n = num_top_strains) %>%
  pull(strain)

combined_rpkmf <- combined_rpkmf %>%
  mutate(strain_grouped = case_when(
    strain == "No Mastadenovirus Detected" ~ "No Mastadenovirus Detected",
    strain %in% top_strains ~ strain,
    TRUE ~ "Other"
  )) %>%
  group_by(sample_ID, Date, strain_grouped) %>%
  summarise(RPKMF = sum(total_RPKMF, na.rm = TRUE), .groups = "drop")

combined_rpkmf <- combined_rpkmf %>%
  group_by(sample_ID, Date) %>%
  mutate(total_sample_RPKMF = sum(RPKMF, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(relative_abundance = ifelse(
    total_sample_RPKMF > 0,
    (RPKMF / total_sample_RPKMF) * 100,
    ifelse(strain_grouped == "No Mastadenovirus Detected", 100, 0)
  ))

# Compute respiratory adenovirus fraction (all except Human adenovirus 41)
sample_sums <- combined_rpkmf %>%
  group_by(sample_ID, Date) %>%
  summarise(
    total_adeno = sum(RPKMF[strain_grouped != "No Mastadenovirus Detected"], 
                      na.rm = TRUE),
    resp_adeno = sum(RPKMF[!strain_grouped %in% 
                             c("No Mastadenovirus Detected", "Human adenovirus 41")], 
                     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(resp_adeno_pct = ifelse(total_adeno > 0, 
                                 (resp_adeno / total_adeno) * 100, 0))

# Aggregate to daily time series
wastewater_timeseries <- sample_sums %>%
  group_by(Date) %>%
  summarise(resp_adeno_pct = mean(resp_adeno_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(Date)

# Load clinical data for adenovirus
clinical_data_adeno <- read_csv("data/public_cases.csv") %>%
  filter(Subject == "V_ADV",
         DateMonday >= as.Date("2022-04-24") & 
           DateMonday <= as.Date("2024-12-01")) %>%
  group_by(Date = DateMonday) %>%
  summarise(cases = n(), .groups = "drop") %>%
  arrange(Date)

# Fuzzy join clinical and wastewater data
joined_data_adeno <- difference_left_join(
  clinical_data_adeno, 
  wastewater_timeseries,
  by = c("Date" = "Date"),
  max_dist = 5,
  distance_col = "date_diff"
) %>%
  group_by(Date.x) %>%
  filter(abs(date_diff) == min(abs(date_diff))) %>%
  ungroup() %>%
  rename(Date = Date.x) %>%
  drop_na() %>%
  arrange(Date)

# Define three time windows for adenovirus analysis
windows <- list(
  A = c("2022-04-25", "2022-08-01"),
  B = c("2022-09-15", "2023-08-03"),
  C = c("2023-10-01", "2024-05-01")
)

window_colors <- list(
  A = "#0fb4d6",
  B = "#91be6d",
  C = "#f5b267"
)

# Perform DTW analysis for each window
dtw_results <- list()
lag_plots <- list()

for(win in names(windows)) {
  win_range <- windows[[win]]
  
  data_win <- joined_data_adeno %>%
    filter(Date >= as.Date(win_range[1]) & Date <= as.Date(win_range[2]))
  
  if(nrow(data_win) < 3) {
    message(paste("Not enough data points in Window", win))
    next
  }
  
  # Compute moving averages for this window
  data_win <- data_win %>%
    arrange(Date) %>%
    mutate(
      smooth_resp = scale(rollapply(
        resp_adeno_pct, 
        width = 2, 
        FUN = function(x) mean(x, na.rm = TRUE),
        fill = NA, 
        align = "right", 
        partial = TRUE
      )),
      smooth_cases = scale(rollapply(
        cases, 
        width = 2, 
        FUN = function(x) mean(x, na.rm = TRUE),
        fill = NA, 
        align = "right", 
        partial = TRUE
      ))
    ) %>%
    drop_na()
  
  # Prepare time series for DTW
  ts_wastewater <- data_win$smooth_resp
  ts_clinical <- data_win$smooth_cases
  
  dtw_alignment <- dtw(
    ts_wastewater, 
    ts_clinical,
    keep = TRUE,
    step.pattern = symmetric2,
    window.type = "sakoechiba",
    window.size = 4
  )
  
  # Compute lags
  lags_dtw <- dtw_alignment$index2 - dtw_alignment$index1
  lag_df <- data.frame(lag = lags_dtw)
  
  # Create complete frequency table
  lag_counts <- lag_df %>%
    count(lag) %>%
    complete(lag = -3:3, fill = list(n = 0))
  
  # Create lag plot
  lag_plot <- ggplot(lag_counts, aes(x = lag, y = n)) +
    geom_col(fill = window_colors[[win]], color = "black", 
             alpha = 0.6, width = 0.9) +
    labs(
      title = paste("Adenovirus - Window", win),
      subtitle = "Negative: Clinical earlier; Positive: WW leads",
      x = "Lag (clinical - wastewater index)",
      y = "Frequency"
    ) +
    scale_x_continuous(breaks = -4:4, limits = c(-5, 5)) +
    scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.25, "cm")
    ) +
    coord_cartesian(clip = "off")
  
  lag_plots[[win]] <- lag_plot
  dtw_results[[win]] <- list(alignment = dtw_alignment, window_data = data_win)
}

# COMBINE ALL LAG PLOTS



combined_lag_panel <- 
  rotavirus_lag_plot + 
  lag_plots[["A"]] + 
  lag_plots[["B"]] + 
  lag_plots[["C"]] +
  plot_layout(ncol = 4) +
  plot_annotation(
    title = "DTW Lag Histograms: Rotavirus & Adenovirus Time Windows",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(combined_lag_panel)

