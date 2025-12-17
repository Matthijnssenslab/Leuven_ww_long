# Primate erythroparvovirus 1 (B19V) Wastewater Surveillance Analysis
# Time series analysis with 4-week moving average

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(zoo)


# Load pre-filtered data
filtered_data <- readRDS("pub_scripts/data/filtered_data_wwlong.rds")

# Filter for Primate erythroparvovirus 1
pep1_data <- filtered_data %>%
  filter(species == "Primate erythroparvovirus 1")

# Sample-level aggregation
sample_agg <- pep1_data %>%
  group_by(sample_ID, Date) %>%
  summarise(RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop")

# Daily aggregation with complete date sequence through end of 2024
daily_data <- sample_agg %>%
  group_by(Date) %>%
  summarise(RPKMF = sum(RPKMF, na.rm = TRUE), .groups = "drop") %>%
  complete(
    Date = seq.Date(as.Date("2022-04-01"), as.Date("2024-12-31"), by = "day"),
    fill = list(RPKMF = 0)
  ) %>%
  arrange(Date) %>%
  mutate(log2_RPKMF = log2(RPKMF + 1))

# Weekly aggregation
weekly_data <- daily_data %>%
  group_by(Week = floor_date(Date, unit = "week", week_start = 1)) %>%
  summarise(weekly_avg = mean(log2_RPKMF, na.rm = TRUE), .groups = "drop") %>%
  arrange(Week)

# Apply 4-week moving average and handle NAs
weekly_data <- weekly_data %>%
  mutate(
    moving_avg = zoo::rollmean(weekly_avg, k = 4, fill = NA, align = "right"),
    moving_avg = ifelse(is.na(moving_avg) | !is.finite(moving_avg), 0, moving_avg)
  )

cat("Total weeks processed:", nrow(weekly_data), "\n")
cat("Date range:", as.character(min(weekly_data$Week)), "to", 
    as.character(max(weekly_data$Week)), "\n")

# Create visualization

pep1_plot <- ggplot(weekly_data, aes(x = Week, y = moving_avg)) +
  geom_area(fill = "#AED6F1", alpha = 0.35) +
  geom_line(color = "#2E86C1", linewidth = 2.2) +
  labs(
    title = "Primate erythroparvovirus 1 in Wastewater Surveillance",
    subtitle = "4-week Moving Average of log2(RPKMF + 1) (2022–2024)",
    x = "Week",
    y = "log2(RPKMF + 1) (Moving Avg)"
  ) +
  scale_x_date(
    date_breaks = "4 months",
    date_labels = "%b\n%Y",
    limits = c(as.Date("2022-04-01"), as.Date("2024-12-31")),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(weekly_data$moving_avg) * 1.05),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 19, color = "#2E86C1", 
                              margin = margin(b = 8)),
    plot.subtitle = element_text(size = 13, color = "black", 
                                 margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(angle = 37, vjust = 0.5, size = 13, 
                               color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 1.5),
    axis.ticks.x = element_line(color = "black", linewidth = 1.5),
    axis.ticks.y = element_line(color = "black", linewidth = 1.5),
    axis.ticks.length.x = unit(0.4, "cm"),
    axis.ticks.length.y = unit(0.4, "cm"),
    plot.margin = margin(18, 18, 18, 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

print(pep1_plot)
