## Figure 2C–D: SARS-CoV-2 clinical incidence vs wastewater signal (no log scale)

# Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  library(scales)
  library(zoo)
})

START_DATE <- as.Date("2022-04-24")
END_DATE   <- as.Date("2024-12-31")

# Wastewater data (pre-filtered dataset, shared with Figure 1)
filtered_data <- readRDS("data/filtered_data_wwlong.rds")

# Weekly SARS-CoV-2 wastewater signal (RPKMF)
ww_sars2_weekly <- filtered_data %>%
  filter(species == "Severe acute respiratory syndrome-related coronavirus") %>%
  mutate(Week = floor_date(Date, "week", week_start = 1)) %>%
  group_by(Week) %>%
  summarise(ww_RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

# Weekly clinical SARS-CoV-2 data (already aggregated)
clin_path <- "data/uzleuven_pathogens_weekly_long.csv"
df_long   <- readr::read_csv(clin_path, show_col_types = FALSE)

clin_sars2_weekly <- df_long %>%
  filter(pathogen %in% c("Coronavirus SARS / COVID-19", 
                         "Coronavirus SARS",
                         "SARS-CoV-2")) %>%
  mutate(Week = floor_date(date, "week", week_start = 1)) %>%
  group_by(Week) %>%
  summarise(clin_cases = sum(amount, na.rm = TRUE), .groups = "drop")

# Merge and complete weekly grid
weeks_all <- tibble(
  Week = seq.Date(
    floor_date(START_DATE, "week", week_start = 1),
    floor_date(END_DATE,   "week", week_start = 1),
    by = "week"
  )
)

aligned <- weeks_all %>%
  left_join(ww_sars2_weekly,  by = "Week") %>%
  left_join(clin_sars2_weekly, by = "Week") %>%
  mutate(
    ww_RPKMF   = replace_na(ww_RPKMF, 0),
    clin_cases = replace_na(clin_cases, 0)
  )

# 4-week moving averages (raw values, not logged)
combined_data <- aligned %>%
  arrange(Week) %>%
  mutate(
    ma_clinical = rollapply(clin_cases, width = 4,
                            FUN = \(x) mean(x, na.rm = TRUE),
                            fill = NA, align = "right"),
    ma_RPKMF    = rollapply(ww_RPKMF, width = 4,
                            FUN = \(x) mean(x, na.rm = TRUE),
                            fill = NA, align = "right")
  )

combined_filtered <- combined_data %>%
  filter(!is.na(ma_clinical) & !is.na(ma_RPKMF))

# Dual-axis scaling (raw values)
range_clin <- range(combined_filtered$ma_clinical, na.rm = TRUE)
range_rpkm <- range(combined_filtered$ma_RPKMF,    na.rm = TRUE)

scale_factor <- diff(range_clin) / diff(range_rpkm)
shift        <- range_clin[1] - range_rpkm[1] * scale_factor

combined_data <- combined_data %>%
  mutate(ma_RPKMF_scaled = ma_RPKMF * scale_factor + shift)

# Panel C — time series (raw values)
SMOOTH_SPAN <- 0.4
SMOOTH_N    <- 600

plot_df <- combined_data %>%
  filter(!is.na(ma_clinical), !is.na(ma_RPKMF_scaled),
         Week >= START_DATE, Week <= END_DATE)

start_date <- min(plot_df$Week)
end_date   <- max(plot_df$Week)
six_months <- seq(from = start_date, to = end_date, by = "6 months")
x_breaks   <- six_months[six_months < end_date]

p_time <- ggplot(plot_df, aes(x = Week)) +
  # ribbons
  geom_smooth(aes(y = ma_clinical, fill = "Clinical cases"),
              method = "loess", span = SMOOTH_SPAN,
              se = TRUE, n = SMOOTH_N, color = NA, alpha = 0.10) +
  geom_smooth(aes(y = ma_RPKMF_scaled, fill = "RPKMF"),
              method = "loess", span = SMOOTH_SPAN,
              se = TRUE, n = SMOOTH_N, color = NA, alpha = 0.10) +
  # dashed raw lines
  geom_line(aes(y = ma_clinical, color = "Clinical cases"),
            linewidth = 1, linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = ma_RPKMF_scaled, color = "RPKMF"),
            linewidth = 1, linetype = "dotdash", alpha = 0.4) +
  # solid smoothed lines
  geom_smooth(aes(y = ma_clinical, color = "Clinical cases"),
              method = "loess", span = SMOOTH_SPAN,
              n = SMOOTH_N, se = FALSE, linewidth = 1.2) +
  geom_smooth(aes(y = ma_RPKMF_scaled, color = "RPKMF"),
              method = "loess", span = SMOOTH_SPAN,
              n = SMOOTH_N, se = FALSE, linewidth = 1.2) +
  scale_y_continuous(
    name = "Clinical cases (4-week MA)",
    sec.axis = sec_axis(~ (. - shift) / scale_factor,
                        name = "RPKMF (4-week MA)")
  ) +
  scale_x_date(
    breaks = x_breaks,
    labels = date_format("%b %Y"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_color_manual(values = c("Clinical cases" = "darkblue",
                                "RPKMF"         = "#ff7f0e")) +
  scale_fill_manual(values  = c("Clinical cases" = "darkblue",
                                "RPKMF"         = "#ff7f0e")) +
  labs(
    title = "SARS-CoV-2: clinical cases vs wastewater signal (raw values, 4-week MA)",
    x = "Week"
  ) +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.title.y.left  = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "#ff7f0e"),
    panel.border       = element_rect(color = "black", fill = NA),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(color = "black")
  )

# Panel D — correlation scatter (raw values)
pearson_ma <- cor.test(
  combined_filtered$ma_RPKMF,
  combined_filtered$ma_clinical,
  method = "pearson"
)

p_scatter <- ggplot(combined_filtered,
                    aes(x = ma_RPKMF, y = ma_clinical)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE,
              color = "darkgreen", linetype = "dashed",
              fill = "darkgreen", alpha = 0.15) +
  labs(
    x = "RPKMF (4-week MA)",
    y = "Clinical cases (4-week MA)"
  ) +
  annotate("text",
           x = min(combined_filtered$ma_RPKMF) + 0.05,
           y = max(combined_filtered$ma_clinical) - 0.05 * max(combined_filtered$ma_clinical),
           hjust = 0,
           label = paste0(
             "Pearson r = ", round(pearson_ma$estimate, 2),
             " (p = ", sprintf("%.4f", pearson_ma$p.value), ")"
           )) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA),
    axis.ticks        = element_line(color = "black")
  )

print(p_time)
print(p_scatter)
