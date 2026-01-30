# Figure 1D: t-SNE by Season


suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
  library(Rtsne)
  library(ggplot2)
  library(ggnewscale)
})

# Load filtered data
filtered_data <- readRDS("data/filtered_data_wwlong.rds")

# Add temporal features
filtered_data <- filtered_data %>%
  mutate(
    Year = year(Date),
    MonthNum = month(Date),
    Season = case_when(
      MonthNum %in% c(12, 1, 2) ~ "Winter",
      MonthNum %in% c(3, 4, 5) ~ "Spring",
      MonthNum %in% c(6, 7, 8) ~ "Summer",
      MonthNum %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Season))

# Create community matrix
community_matrix <- filtered_data %>%
  group_by(sample_ID, species) %>%
  summarise(RPKMF = log10(sum(RPKMF, na.rm = TRUE) + 1), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = RPKMF, values_fill = 0) %>%
  column_to_rownames("sample_ID") %>%
  as.matrix()

if (nrow(community_matrix) < 5) {
  stop("Not enough samples after filtering to compute t-SNE.")
}

# PCA → t-SNE
set.seed(6)
num_samples <- nrow(community_matrix)
perp <- max(5, min(30, floor(num_samples / 3)))

pca_results <- prcomp(community_matrix, center = TRUE, scale. = FALSE)
num_pcs <- min(10, ncol(pca_results$x))

if (num_pcs < 2) {
  stop("Not enough components for t-SNE. Check input data.")
}

tsne_results <- Rtsne(
  pca_results$x[, 1:num_pcs, drop = FALSE],
  perplexity = perp,
  check_duplicates = FALSE,
  pca = FALSE,
  verbose = FALSE
)

# Build plotting dataframe
tsne_df <- data.frame(
  tSNE1 = tsne_results$Y[, 1],
  tSNE2 = tsne_results$Y[, 2],
  sample_ID = rownames(community_matrix)
) %>%
  left_join(distinct(filtered_data, sample_ID, Date, Season, Year), by = "sample_ID") %>%
  mutate(
    Year = factor(Year),
    YearProgress = (yday(Date) - 1) / (if_else(leap_year(Date), 366, 365) - 1)
  )

# Season centroids
season_centroids <- tsne_df %>%
  group_by(Season) %>%
  summarise(c1 = mean(tSNE1, na.rm = TRUE), c2 = mean(tSNE2, na.rm = TRUE), .groups = "drop")

# Season colors
season_colors <- c(
  "Winter" = "#4C6EF5",
  "Spring" = "#2FB344",
  "Summer" = "#F59F00",
  "Fall"   = "#D9480F"
)

# Create plot
p <- ggplot(tsne_df, aes(tSNE1, tSNE2)) +
  stat_ellipse(aes(fill = Season), type = "t", level = 0.8,
               geom = "polygon", alpha = 0.12, color = NA) +
  scale_fill_manual(values = season_colors, drop = TRUE, guide = "none") +
  ggnewscale::new_scale_color() +
  geom_point(aes(color = YearProgress, shape = Year), size = 2.6, alpha = 0.95) +
  scale_color_gradientn(
    colors = c("#4C6EF5", "#2FB344", "#F59F00", "#D9480F", "#4C6EF5"),
    values = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("Jan", "Apr", "Jul", "Oct", "Dec"),
    name = "Month (Jan→Dec)",
    guide = guide_colorbar(
      ticks = TRUE,
      draw.ulim = TRUE,
      draw.llim = TRUE,
      barwidth = unit(0.025, "npc"),
      barheight = unit(0.45, "npc")
    )
  ) +
  geom_point(
    data = season_centroids,
    aes(x = c1, y = c2, fill = Season),
    shape = 21, color = "black", stroke = 0.3, size = 4, alpha = 0.8,
    inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    legend.position = "right",
    legend.key.height = unit(0.5, "lines"),
    legend.key.width  = unit(1.0, "lines"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  ) +
  labs(
    title = "t-SNE: Seasons (ellipses/centroids), Color=Month progression, Shape=Year",
    x = "t-SNE 1",
    y = "t-SNE 2",
    shape = "Year"
  )

print(p)

ggsave("tsne_all_seasons.png", plot = p, width = 9, height = 7, dpi = 300)
cat("Saved: tsne_all_seasons.png\n")

# dbRDA: cyclic month effect
pred <- filtered_data %>%
  distinct(sample_ID, Date, Year) %>%
  mutate(
    days_in_year = if_else(leap_year(Date), 366, 365),
    year_frac = (yday(Date) - 0.5) / days_in_year,
    theta = 2*pi*year_frac,
    sin_month = sin(theta),
    cos_month = cos(theta),
    Year = factor(Year)
  ) %>%
  filter(sample_ID %in% rownames(community_matrix)) %>%
  arrange(match(sample_ID, rownames(community_matrix)))

cap <- capscale(
  community_matrix ~ sin_month + cos_month + Condition(Year),
  data = pred,
  distance = "bray"
)

anova_overall <- anova.cca(cap, permutations = 999)
anova_terms   <- anova.cca(cap, by = "terms", permutations = 999)
r2_adj        <- RsquareAdj(cap)

print(anova_overall)
print(anova_terms)
print(r2_adj)

# Save summary
summary_lines <- c(
  "dbRDA with cyclic month predictors (sin/cos), Condition on Year",
  "---------------------------------------------------------------",
  capture.output(anova_overall),
  "",
  capture.output(anova_terms),
  "",
  paste("Adjusted R2:", round(r2_adj$adj.r.squared, 4))
)
writeLines(summary_lines, "dbrda_month_results.txt")
cat("Saved: dbrda_month_results.txt\n")
