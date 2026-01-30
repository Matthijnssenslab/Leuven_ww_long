# Figure 1C: Rarefaction Curves - Species vs Strain


suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(micropan)
  library(wesanderson)
  library(tibble)
})

utils::globalVariables(c("taxon", "value"))

# Load filtered data
filtered_data <- readRDS("data/filtered_data_wwlong.rds") %>%
  mutate(value_for_matrix = if ("RPKMF" %in% names(.)) .data[["RPKMF"]] else .data[["reads_aligned"]])

# Compute rarefaction
compute_rarefaction <- function(long_dt, taxon_col, label) {
  taxon <- NULL; value <- NULL
  mat <- long_dt %>%
    group_by(.data[["sample_ID"]], taxon = .data[[taxon_col]]) %>%
    summarise(value = mean(.data[["value_for_matrix"]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = value, values_fill = 0) %>%
    column_to_rownames(var = "sample_ID") %>%
    as.matrix()
  if (ncol(mat) > 0) mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  rr <- rarefaction(mat, n.perm = 100)
  rr$avg   <- rowMeans(rr[, -1, drop = FALSE], na.rm = TRUE)
  rr$min_g <- apply(rr[, -1, drop = FALSE], 1, min, na.rm = TRUE)
  rr$max_g <- apply(rr[, -1, drop = FALSE], 1, max, na.rm = TRUE)
  rr$type <- label
  rr
}

rare_species <- compute_rarefaction(filtered_data, "species", "Species")
rare_strain  <- compute_rarefaction(filtered_data, "strain",  "Strain")
rare_both <- bind_rows(rare_species, rare_strain)

# Create plot
pal <- wes_palette("Darjeeling1", type = "continuous")
col_map <- c("Species" = pal[2], "Strain" = pal[1])

p <- ggplot(rare_both, aes(x = Genome, y = avg, color = type, fill = type)) +
  geom_ribbon(aes(ymin = min_g, ymax = max_g, group = type), alpha = 0.25, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  labs(
    title = "Rarefaction Curves: Species vs Strain",
    x = "Number of Samples Analyzed",
    y = "Unique Taxa Detected",
    color = "Level",
    fill = "Level"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

print(p)

# Print summary statistics
num_species <- filtered_data %>% distinct(species) %>% nrow()
num_strains <- filtered_data %>% distinct(strain) %>% nrow()
num_samples <- filtered_data %>% distinct(sample_ID) %>% nrow()
cat("Unique species:", num_species, "\n")
cat("Unique strains:", num_strains, "\n")
cat("Unique samples:", num_samples, "\n")

