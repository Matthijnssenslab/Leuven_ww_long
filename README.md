# Hybrid-capture-enabled longitudinal metagenomics allows strain-resolved human-associated virus surveillance in wastewater

This repository contains processed datasets and R scripts to reproduce the analyses and figures from the study. It also documents how the raw reads were processed so the workflow can be reproduced end-to-end.

## Dataset (raw reads)

BioProject: PRJNA1391037 (SRA). Human-associated reads were removed before deposition.

Raw reads can be downloaded using the SRA Toolkit (e.g., `prefetch` / `fasterq-dump`). See:
https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

Processed inputs for figure reproduction are included in `data/` so you can regenerate figures without re-running the raw-read pipeline.

## Metagenomics analyses (EsViritu v0.2.3)

Raw reads were processed with EsViritu v0.2.3 using the Virus Pathogen Database v2.0.2 (GenBank content through Nov 2022; Zenodo 7876309) as the reference.
EsViritu links:
- GitHub: https://github.com/cmmr/EsViritu/
- Paper (Tisza et al., Nat Commun 2023): https://www.nature.com/articles/s41467-023-42064-1
If you use EsViritu, please cite the related paper.

Workflow summary:
- Quality filtering and adapter trimming with fastp, with deduplication enabled.
- Reference mapping at >=90% nucleotide identity and >=90% read coverage.
- Consensus sequence generation with samtools; near-duplicate consensus sequences removed at >95% similarity.
- The pipeline reports all findings (we have removed internal reporting cutoff of EsViritu); downstream filtering was done in RStudio with cutoffs mentioned below.

These steps produce the processed datasets in `data/` that are used by the figure scripts.

## Pre-filtering in R (create_filtered_data.R)

The pre-filtering step is implemented in `scripts/create_filtered_data.R`. It removes known wet-lab contaminants and low-confidence hits to improve certainty of detections and downstream genomic assignments. We also exclude a short-amplicon contamination in Human mastadenovirus C (<2 kb). Two filtered datasets are produced:
- Standard: >=1000 bp coverage or >=50% completeness, and >10 reads.
- Relaxed: >=500 bp coverage or >=50% completeness, and >10 reads, used for analyses where low-read targets would otherwise be lost; after contaminant removal and the read-count filter, the chance of spurious low-coverage hits is reduced.

## Clinical datasets

Clinical datasets used in the study are included in this repository for reproducibility. These data sources are cited in the manuscript and are available under `data/` (e.g., `public_cases.csv`, `rva_weekly_cases.csv`, `uzleuven_pathogens_weekly_long.csv`).

## Rotavirus A genotyping

Rotavirus A segment genotypes were assigned by competitive mapping against RCWG genotype reference sets:
- Reference FASTA/TSV files per segment are stored in `data/database_for_RVA_genotyping/`.
- The competitive mapping workflow used `influenza_a_serotype` v0.1.5, replacing its influenza references with the RCWG segment references.
- For each segment, reads are aligned to all genotype references and the best-supported genotype is assigned; ambiguous reads are marked as ambiguous and are not shown in figures.

influenza_a_serotype links and related paper:
- GitHub: https://github.com/mtisza1/influenza_a_serotype
- Related paper (Tisza et al., N Engl J Med 2024): https://www.nejm.org/doi/full/10.1056/NEJMc2405937
If you use influenza_a_serotype, please cite the related paper.

## Repository structure

```
.
├── data/                       # Processed datasets (RDS and CSV)
│   ├── filtered_data_wwlong.rds
│   ├── filtered_data_wwlong_500bp.rds
│   ├── processed_clinical_rotavirus.rds
│   ├── processed_wastewater_rotavirus.rds
│   ├── public_cases.csv
│   └── ...
├── scripts/                    # R scripts for generating figures
│   ├── Figure_1A.R             # Percentage of aligned reads
│   ├── Figure_1B.R             # Sequences per coverage category
│   ├── Figure_3ABCDEFGHI.R     # Rotavirus genotype analysis
│   └── ...
├── figures_pdf/                # Output figures in PDF format
└── README.md                   # Project documentation
```

## Prerequisites

System requirements:
- Operating system: macOS (analysis was run on macOS)
- R version: >= 4.0.0
- No special hardware required.

R dependencies:
```r
install.packages(c(
  "tidyverse",
  "lubridate",
  "ggplot2",
  "patchwork",
  "RColorBrewer",
  "gridExtra",
  "scales",
  "grid",
  "ggpubr",
  "ggrepel",
  "viridis",
  "zoo"
))
```

## Usage (figure reproduction)

1. Clone the repository:
   ```bash
   git clone https://github.com/Matthijnssenslab/Leuven_ww_long.git
   cd Leuven_ww_long
   ```

2. Set the working directory in R/RStudio to the repository root:
   ```r
   setwd("/path/to/Leuven_ww_long")
   ```

3. Run scripts named after the figures they generate, for example:
   ```r
   source("scripts/Figure_1A.R")
   ```

Note: Some figure panels were generated separately and assembled in Adobe Illustrator 2024 without editing data-associated elements. Virus taxonomy labels were updated to current ICTV naming. Run time is <1 min for each figure.

## Data description

- `filtered_data_wwlong.rds`: Main dataset containing filtered sequencing reads.
- `processed_clinical_rotavirus.rds`: Anonymized clinical surveillance data for rotavirus comparisons.
- `processed_wastewater_rotavirus.rds`: Aggregated wastewater data for rotavirus genotyping.

## License

MIT License. See `LICENSE`.

## Citation

If you use this code or data, please cite:
> Mustafa Karatas, Mandy Bloemen, Jill Swinnen, Lila Close, Marc Van Ranst, Elke Wollants, Jelle Matthijnssens. "Hybrid-capture-enabled longitudinal metagenomics allows strain-resolved human-associated virus surveillance in wastewater". SSRN Preprint. DOI: https://doi.org/10.2139/ssrn.6047769
