## Figure 2A: Rotavirus clinical vs wastewater genomics (weekly)

# Load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyquant)
library(lubridate)
library(scales)
library(viridis)
library(car)
library(ggpubr)
library(broom)
library(stringr)
library(zoo)
library(grid)

# Wastewater data (pre-filtered, sample-level)
# Uses data/filtered_data_wwlong.rds (same filtering as Figure 1)
filtered_data <- readRDS("data/filtered_data_wwlong.rds") %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= as.Date("2022-04-24") & Date <= as.Date("2024-12-31"))

# Weekly aggregation of Rotavirus RPKMF
weekly_rpkmf <- filtered_data %>%
  filter(species == "Rotavirus A") %>%
  mutate(Week = floor_date(Date, "week", week_start = 1)) %>%
  group_by(Week) %>%
  summarise(mean_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop") %>%
  complete(
    Week = seq.Date(min(Week), max(Week), by = "week"),
    fill = list(mean_RPKMF = NA_real_)
  )

# Clinical data (already aggregated weekly)
# Uses data/rva_weekly_cases.csv with:
#   DateMonday (week start, YYYY-MM-DD)
#   clinical_cases (weekly count)
weekly_clinical <- readr::read_csv("data/rva_weekly_cases.csv", show_col_types = FALSE) %>%
  mutate(DateMonday = as.Date(DateMonday)) %>%
  filter(DateMonday >= as.Date("2022-04-01") & DateMonday <= as.Date("2024-12-31")) %>%
  complete(
    DateMonday = seq.Date(min(DateMonday), max(DateMonday), by = "week"),
    fill = list(clinical_cases = 0L)
  )

# Combine weekly clinical and wastewater data
combined_data <- inner_join(
  weekly_clinical,
  weekly_rpkmf,
  by = c("DateMonday" = "Week")
) %>%
  arrange(DateMonday)

# Log-transform with +1 guards (handle zeros) for plotting/MA
combined_data <- combined_data %>%
  mutate(
    log_clinical_cases = log2(clinical_cases + 1),
    log_mean_RPKMF = log10(mean_RPKMF + 1)
  )

# 4-week moving averages (on log-transformed metrics)
combined_data <- combined_data %>%
  arrange(DateMonday) %>%
  mutate(
    ma_clinical = zoo::rollapply(log_clinical_cases, width = 4,
                                 FUN = function(x) mean(x, na.rm = TRUE),
                                 fill = NA, align = "right", partial = FALSE),
    ma_RPKMF = zoo::rollapply(log_mean_RPKMF, width = 4,
                              FUN = function(x) mean(x, na.rm = TRUE),
                              fill = NA, align = "right", partial = FALSE)
  )

combined_data_filtered <- combined_data %>%
  filter(!is.na(ma_clinical) & !is.na(ma_RPKMF))

# Scale RPKMF MA to left axis
range_clin <- range(combined_data$ma_clinical, na.rm = TRUE)
range_rpkm <- range(combined_data$ma_RPKMF, na.rm = TRUE)
scale_factor <- diff(range_clin) / diff(range_rpkm)
shift <- range_clin[1] - range_rpkm[1] * scale_factor

combined_data <- combined_data %>%
  mutate(ma_RPKMF_scaled = ma_RPKMF * scale_factor + shift)


# Smoother settings
SMOOTH_SPAN <- 0.40
SMOOTH_N    <- 600

# Time series plot with dual y-axes

plot_df <- combined_data %>%
  dplyr::filter(!is.na(ma_clinical), !is.na(ma_RPKMF_scaled))

# X-axis breaks: start, then every 6 months
start_date <- min(plot_df$DateMonday, na.rm = TRUE)
end_date   <- max(plot_df$DateMonday, na.rm = TRUE)

sixmo_seq  <- seq(from = start_date, to = end_date, by = "6 months")
x_breaks   <- sort(unique(c(start_date, sixmo_seq)))
x_breaks   <- x_breaks[x_breaks < end_date]  # drop end date if present

p_time_series <- ggplot(plot_df, aes(x = DateMonday)) +
  geom_smooth(aes(y = ma_clinical, fill = "Clinical Cases"),
              method = "loess", span = SMOOTH_SPAN, se = TRUE, n = SMOOTH_N,
              color = NA, alpha = 0.10) +
  geom_smooth(aes(y = ma_RPKMF_scaled, fill = "RPKMF"),
              method = "loess", span = SMOOTH_SPAN, se = TRUE, n = SMOOTH_N,
              color = NA, alpha = 0.10) +
  
  geom_line(aes(y = ma_clinical, color = "Clinical Cases"),
            linewidth = 1.0, linetype = "dashed", alpha = 0.40, na.rm = TRUE) +
  
  geom_line(aes(y = ma_RPKMF_scaled, color = "RPKMF"),
            linewidth = 1.0, linetype = "dotdash", alpha = 0.40, na.rm = TRUE) +
  
  geom_smooth(aes(y = ma_clinical, color = "Clinical Cases"),
              method = "loess", span = SMOOTH_SPAN, se = FALSE, n = SMOOTH_N,
              linewidth = 1.2) +
  
  geom_smooth(aes(y = ma_RPKMF_scaled, color = "RPKMF"),
              method = "loess", span = SMOOTH_SPAN, se = FALSE, n = SMOOTH_N,
              linewidth = 1.2) +
  
  scale_y_continuous(
    name = "log2(Clinical Cases + 1)",
    sec.axis = sec_axis(~ (. - shift) / scale_factor, name = "log10(RPKMF + 1)")
  ) +
  scale_x_date(
    breaks = x_breaks,
    labels = scales::date_format("%b %Y"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_color_manual(values = c("Clinical Cases" = "darkblue",
                                "RPKMF" = "#ff7f0e")) +
  scale_fill_manual(values  = c("Clinical Cases" = "darkblue",
                                "RPKMF" = "#ff7f0e")) +
  labs(
    title = "Rotavirus: log2(Clinical Cases) vs. log10(RPKMF) (4-week Moving Average)",
    x = "Date", color = "Metric", fill = NULL
  ) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    axis.title.y.left  = element_text(color = "darkblue", size = 12),
    axis.title.y.right = element_text(color = "#ff7f0e", size = 12),
    panel.border       = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(color = "black")
  )

print(p_time_series)


# Correlation scatter on 4-week moving averages
pearson_ma <- cor.test(combined_data_filtered$ma_RPKMF,
                       combined_data_filtered$ma_clinical,
                       method = "pearson", use = "complete.obs")

p_scatter <- ggplot(combined_data_filtered,
                    aes(x = ma_RPKMF, y = ma_clinical)) +
  geom_point(shape = 16, stroke = 0, alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, n = 500,
              color = "darkgreen", linetype = "dashed",
              fill = "darkgreen", alpha = 0.15) +
  labs(title = "Correlation on 4-week Moving Averages",
       x = "log10(RPKMF + 1) (4-week MA)",
       y = "log2(Clinical Cases + 1) (4-week MA)") +
  annotate("text",
           x = min(combined_data_filtered$ma_RPKMF, na.rm = TRUE) + 0.05,
           y = max(combined_data_filtered$ma_clinical, na.rm = TRUE) - 0.1,
           hjust = 0,
           label = paste0("Pearson r = ", round(pearson_ma$estimate, 2),
                          " (p = ", sprintf("%.4f", pearson_ma$p.value), ")")) +
  theme_minimal() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p_scatter)

