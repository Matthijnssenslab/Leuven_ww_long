# Create Filtered Datasets for wwlong Scripts
# Run this script first to create filtered_data_wwlong.rds files
# Rationale: remove known contaminants and low-confidence hits to improve detection certainty.
# A relaxed 500 bp threshold is used for some analyses to retain low-read targets.

setwd("path/to/your/tsv/from/EsViritu")

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

metadata <- readr::read_csv(
  "csv/WW_long_metadata.csv",
  col_types = cols(Date = col_character())) %>%
  mutate(Date = suppressWarnings(as.Date(Date, format = "%Y-%m-%d"))) %>%
  filter(!is.na(Date))

data_tsv <- readr::read_tsv("wastewater_long_full.detected_virus.combined.tax.tsv", col_types = cols())

# Filter for samples in metadata (all samples are Belgium)
belgium_sample_ids <- metadata$sample_ID

merged_data <- data_tsv %>%
  filter(sample_ID %in% belgium_sample_ids) %>%
  left_join(metadata, by = "sample_ID")

# Standard filtered dataset (1000bp threshold)
filtered_data <- merged_data %>%
  # Remove known wet-lab contaminants observed previously:
  # Parvovirus NIH-CQV, Chikungunya virus, unclassified Culex Bastrovirus-like virus,
  # unclassified Mus musculus mobilized endogenous polytropic provirus.
  filter(
    !str_detect(strain, "Parvovirus NIH-CQV"),
    !str_detect(strain, "Chikungunya virus"),
    !str_detect(strain, "unclassified Culex Bastrovirus-like virus strain"),
    !str_detect(strain, "unclassified Mus musculus mobilized endogenous polytropic provirus strain"),
    # Amplicon contamination: remove Human mastadenovirus C <2 kb.
    !(species == "Human mastadenovirus C" & covered_bases <= 2000)
  ) %>%
  # Remove low-certainty hits: >=1000 bp OR >=50% completeness and >10 reads.
  filter(
    (covered_bases >= 1000 | (covered_bases / reference_length) * 100 >= 50) &
      reads_aligned > 10
  ) %>%
  # Remove additional known contaminants (DCCV-4, DCCV-13, HERV-K) and low-depth samples.
  filter(
    !(strain %in% c("Circovirus-like genome DCCV-4", "Human endogenous retrovirus K", "Circovirus-like genome DCCV-13")),
    total_filtered_reads_in_sample > 1e6
  ) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= as.Date("2022-04-23") & Date <= as.Date("2024-12-31"))

# Relaxed filtered dataset (500bp threshold)
filtered_data_500bp <- merged_data %>%
  # Remove known wet-lab contaminants observed previously:
  # Parvovirus NIH-CQV, Chikungunya virus, unclassified Culex Bastrovirus-like virus,
  # unclassified Mus musculus mobilized endogenous polytropic provirus.
  filter(
    !str_detect(strain, "Parvovirus NIH-CQV"),
    !str_detect(strain, "Chikungunya virus"),
    !str_detect(strain, "unclassified Culex Bastrovirus-like virus strain"),
    !str_detect(strain, "unclassified Mus musculus mobilized endogenous polytropic provirus strain"),
    # Amplicon contamination: remove Human mastadenovirus C <2 kb.
    !(species == "Human mastadenovirus C" & covered_bases <= 2000)
  ) %>%
  # Remove low-certainty hits: >=500 bp OR >=50% completeness and >10 reads.
  filter(
    (covered_bases >= 500 | (covered_bases / reference_length) * 100 >= 50) &
      reads_aligned > 10
  ) %>%
  # Remove additional known contaminants (DCCV-4, DCCV-13, HERV-K) and low-depth samples.
  filter(
    !(strain %in% c("Circovirus-like genome DCCV-4", "Human endogenous retrovirus K", "Circovirus-like genome DCCV-13")),
    total_filtered_reads_in_sample > 1e6
  ) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= as.Date("2022-04-23") & Date <= as.Date("2024-12-31"))

# Save filtered data
if (!dir.exists("pub_scripts/data")) {
  dir.create("pub_scripts/data", recursive = TRUE)
}

saveRDS(filtered_data, "pub_scripts/data/filtered_data_wwlong.rds")
saveRDS(filtered_data_500bp, "pub_scripts/data/filtered_data_wwlong_500bp.rds")

cat("Filtered datasets saved:\n")
cat("  - Standard (1000bp): pub_scripts/data/filtered_data_wwlong.rds\n")
cat("  - Relaxed (500bp):   pub_scripts/data/filtered_data_wwlong_500bp.rds\n")
