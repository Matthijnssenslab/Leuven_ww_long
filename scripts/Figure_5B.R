# Rotavirus "Other" Genotypes, 5B: Co-occurrence Network Analysis
# Creates chord diagram showing genotype co-occurrence patterns with optimal ordering

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(circlize)
library(seriation)


# Load processed rotavirus data
wastewater_data <- readRDS("processed_wastewater_rotavirus.rds")

# Define host origin mapping for genotypes (ONLY these will be included for co-occurence analyses for simplicity.) 
#If you include DS-1-like genotypes or Wa-like genotypes,they will be connected to everything since they are detected throughout the study.
host_origin <- tribble(
  ~genotype, ~host_origin,
  "P5", "Bovine",
  "P11", "Bovine",
  "P14", "Bovine",
  "P9", "Feline/Canine",
  "I3", "Feline/Canine",
  "I20", "Rodents",
  "A3", "Bovine/Feline/Canine",
  "A11", "Bovine",
  "A13", "Bovine",
  "N3", "Bovine",
  "N4", "Avian",
  "T4", "Avian",
  "E4", "Avian",
  "H4", "Avian",
  "R4", "Avian",
  "C4", "Avian",
  "I4", "Avian",
  "T3", "Feline/Canine",
  "T6", "Bovine",
  "T14", "Rodents",
  "E3", "Feline/Canine",
  "H3", "Bovine/Feline/Canine",
  "H6", "Feline/Canine",
  "H13", "Rodents",
  "E12", "Bovine",
  "R3", "Feline/Canine",
  "C3", "Feline/Canine",
  "M3", "Feline/Canine",
  "R11", "Rodents",
  "C11", "Rodents",
  "M10", "Rodents",
  "G6", "Bovine"
)

# Get list of valid genotypes
valid_genotypes <- host_origin$genotype

# Define host-specific colors
host_colors <- c(
  "Bovine" = "#E63946",
  "Feline/Canine" = "#457B9D",
  "Rodents" = "black",
  "Avian" = "darkgreen",
  "Bovine/Feline/Canine" = "#F77F00",
  "Other" = "#CCCCCC"
)

# Prepare co-occurrence data - FILTER for valid genotypes only

chord_data_cooccur <- wastewater_data %>%
  filter(reads > 5) %>%
  filter(geno_group == "other") %>%
  filter(genotype %in% valid_genotypes) %>%  # FILTER: Only include mapped genotypes
  select(sample_ID, segment, genotype) %>%
  distinct() %>%
  filter(segment %in% c("vp7", "vp4", "nsp1", "nsp3", "nsp4", "nsp5"))

cat("Genotypes included:", length(unique(chord_data_cooccur$genotype)), "\n")

# Create all pairwise co-occurrence combinations
pair_list_cooccur <- list()

for (sample in unique(chord_data_cooccur$sample_ID)) {
  genotypes <- chord_data_cooccur %>%
    filter(sample_ID == sample) %>%
    select(segment, genotype) %>%
    distinct() %>%
    mutate(genotype_label = paste0(segment, "_", genotype)) %>%
    pull(genotype_label)
  
  if (length(genotypes) > 1) {
    pairs <- combn(sort(unique(genotypes)), 2, simplify = FALSE)
    for (pair in pairs) {
      pair_list_cooccur[[length(pair_list_cooccur) + 1]] <- data.frame(
        from = pair[1],
        to = pair[2],
        weight = 1,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Count co-occurrences
chord_long_cooccur <- bind_rows(pair_list_cooccur) %>%
  group_by(from, to) %>%
  summarise(weight = n(), .groups = "drop") %>%
  arrange(desc(weight))

# Filter by minimum weight threshold
min_weight_cooccur <- 10
chord_long_cooccur_filtered <- chord_long_cooccur %>%
  filter(weight >= min_weight_cooccur)

cat("Total genotype pairs (≥", min_weight_cooccur, "co-occurrences):", 
    nrow(chord_long_cooccur_filtered), "\n")



# Get all unique genotypes
all_genotypes_cooccur <- unique(c(
  chord_long_cooccur_filtered$from, 
  chord_long_cooccur_filtered$to
))

cat("Unique genotypes in network:", length(all_genotypes_cooccur), "\n")

# Build adjacency matrix for seriation

adj_matrix <- matrix(
  0,
  nrow = length(all_genotypes_cooccur),
  ncol = length(all_genotypes_cooccur),
  dimnames = list(all_genotypes_cooccur, all_genotypes_cooccur)
)

for (i in 1:nrow(chord_long_cooccur_filtered)) {
  g1 <- as.character(chord_long_cooccur_filtered$from[i])
  g2 <- as.character(chord_long_cooccur_filtered$to[i])
  wt <- chord_long_cooccur_filtered$weight[i]
  adj_matrix[g1, g2] <- wt
  adj_matrix[g2, g1] <- wt
}

# Optimal leaf ordering for best arrangement
dist_matrix <- as.dist(max(adj_matrix) - adj_matrix)
order_obj <- seriate(dist_matrix, method = "OLO")
best_order <- get_order(order_obj)
ordered_genotypes_optimal <- all_genotypes_cooccur[best_order]

# Create genotype mapping with host colors and display labels
genotype_mapping_cooccur <- data.frame(
  genotype_label = ordered_genotypes_optimal
) %>%
  separate(genotype_label, into = c("segment", "genotype"), 
           sep = "_", remove = FALSE) %>%
  left_join(host_origin %>% mutate(genotype = as.character(genotype)), 
            by = "genotype") %>%
  mutate(
    color = case_when(
      host_origin == "Bovine" ~ host_colors["Bovine"],
      host_origin == "Feline/Canine" ~ host_colors["Feline/Canine"],
      host_origin == "Rodents" ~ host_colors["Rodents"],
      host_origin == "Avian" ~ host_colors["Avian"],
      host_origin == "Bovine/Feline/Canine" ~ host_colors["Bovine/Feline/Canine"],
      TRUE ~ host_colors["Other"]
    ),
    display_label = genotype
  )

colors_cooccur <- setNames(
  genotype_mapping_cooccur$color, 
  genotype_mapping_cooccur$genotype_label
)

display_labels <- setNames(
  genotype_mapping_cooccur$display_label, 
  genotype_mapping_cooccur$genotype_label
)

# Reorder edge list according to optimal ordering
chord_long_cooccur_filtered <- chord_long_cooccur_filtered %>%
  mutate(
    from = factor(from, levels = ordered_genotypes_optimal),
    to = factor(to, levels = ordered_genotypes_optimal)
  ) %>%
  arrange(from, to)

# Create chord diagram
cat("Creating chord diagram...\n")

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 3,
  track.margin = c(-0.1, 0.1),
  points.overflow.warning = FALSE
)

par(mar = rep(0, 4))

chordDiagram(
  x = chord_long_cooccur_filtered,
  grid.col = colors_cooccur,
  transparency = 0.2,
  directional = 0,
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.05),
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  small.gap = 1
)

# Add genotype labels
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    lab <- display_labels[sector.index]
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 2.5,
      labels = lab,
      facing = "bending",
      niceFacing = TRUE,
      cex = 0.85,
      font = 2
    )
  }
)

# Add legend
legend(
  "topright",
  legend = c("Bovine", "Feline/Canine", "Rodents", "Avian", 
             "Bovine/Feline/Canine", "Other"),
  fill = c("#E63946", "#457B9D", "black", "darkgreen", "#F77F00", "#CCCCCC"),
  border = "black",
  bty = "n",
  cex = 0.85,
  title = "Host Origin",
  title.font = 2
)

#supplementary for the occurence of unusual genotypes

# ============================================================================
# SUPPLEMENTARY BARPLOT: All Relevant Genotype Frequencies (Ordered by Segment)
# ============================================================================


# Get all samples with valid genotypes (5 read minimum, no co-occurrence filter)
all_genotype_data <- wastewater_data %>%
  filter(reads > 5) %>%
  filter(geno_group == "other") %>%
  filter(genotype %in% valid_genotypes) %>%
  select(sample_ID, segment, genotype) %>%
  distinct()

# Define segment order
segment_order <- c("vp7", "vp4", "vp6", "vp1", "vp2", "vp3", 
                   "nsp1", "nsp2", "nsp3", "nsp4", "nsp5")

# Count number of samples per genotype-segment combination
genotype_frequency_all <- all_genotype_data %>%
  group_by(segment, genotype) %>%
  summarise(n_samples = n_distinct(sample_ID), .groups = "drop") %>%
  left_join(host_origin %>% mutate(genotype = as.character(genotype)), 
            by = "genotype") %>%
  mutate(
    host_origin = replace_na(host_origin, "Other"),
    color = case_when(
      host_origin == "Bovine" ~ host_colors["Bovine"],
      host_origin == "Feline/Canine" ~ host_colors["Feline/Canine"],
      host_origin == "Rodents" ~ host_colors["Rodents"],
      host_origin == "Avian" ~ host_colors["Avian"],
      host_origin == "Bovine/Feline/Canine" ~ host_colors["Bovine/Feline/Canine"],
      TRUE ~ host_colors["Other"]
    ),
    # Factor segment for ordering
    segment = factor(segment, levels = segment_order),
    # Create display label with segment info
    genotype_label = paste0(genotype, " (", toupper(segment), ")")
  ) %>%
  arrange(segment, desc(n_samples)) %>%
  # Create ordering variable: order by segment first, then by n_samples within segment
  mutate(plot_order = row_number(),
         genotype_label = factor(genotype_label, levels = genotype_label))

cat("Total genotype-segment combinations:", nrow(genotype_frequency_all), "\n")

# Create barplot with segment information
fig_frequency_all <- ggplot(genotype_frequency_all, 
                            aes(x = genotype_label, 
                                y = n_samples, 
                                fill = host_origin)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_samples), hjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = host_colors, name = "Host Origin") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Frequency of 'Other' Rotavirus Genotypes Across Samples",
    subtitle = "Number of samples where genotype was detected (min 5 reads, ordered by segment then frequency)",
    x = "Genotype (Segment)",
    y = "Number of Samples"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50", 
                                 margin = margin(b = 12)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title = element_text(face = "bold", size = 11),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(3, "pt"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black", 
                                     linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(fig_frequency_all)

# Summary statistics by segment

cat("  Total genotype-segment combinations:", nrow(genotype_frequency_all), "\n")
cat("  Mean samples per genotype:", 
    round(mean(genotype_frequency_all$n_samples), 1), "\n")
cat("  Median samples:           ", 
    median(genotype_frequency_all$n_samples), "\n")
cat("\n  Genotypes by segment (ordered):\n")

segment_summary <- genotype_frequency_all %>%
  group_by(segment) %>%
  summarise(
    n_genotypes = n(), 
    total_detections = sum(n_samples),
    most_common = genotype[which.max(n_samples)],
    max_detections = max(n_samples),
    .groups = "drop"
  )

for (i in 1:nrow(segment_summary)) {
  cat("    ", toupper(as.character(segment_summary$segment[i])), ": ", 
      segment_summary$n_genotypes[i], " genotypes, ",
      segment_summary$total_detections[i], " total detections",
      " (most common: ", segment_summary$most_common[i], " - ",
      segment_summary$max_detections[i], " samples)\n", sep = "")
}

