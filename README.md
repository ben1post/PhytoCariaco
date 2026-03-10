# PhytoCariaco: Phytoplankton Community Dynamics in the Cariaco Basin

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18942210.svg)](https://doi.org/10.5281/zenodo.18942210)

This repository contains the R code used for the data analysis and figure generation presented in the manuscript:

**"Reversible Regime Change: Climate-driven phytoplankton community shifts in the Cariaco Basin, Venezuela"**  
*Under review at JGR Biogeosciences*

Authored by Post, Acevedo-Trejos, Chakraborty, Barton, and Merico.  
All code was written by Post.

---

## Overview

The analysis investigates long-term changes in phytoplankton community structure in the Cariaco Basin using time series analysis and statistics, hierarchical clustering, ordination methods, and an analysis of potential drivers of change using the gradient forest algorithm.

---

## Repository Structure

```
PhytoCariaco/
├── PhytoCariaco.Rproj                  # RStudio project file
├── data/                               # Raw data, processing pipeline, and cleaned outputs
└── plots/                              # Analysis and figure scripts
```

### `data/`

Raw data sources and the full processing pipeline to produce cleaned datasets.

**Raw data sources:**

| Folder | Contents |
|---|---|
| `BCO-DMO/` | Raw CARIACO time series data: `phytoplankton.csv`, `niskin.csv`, `ctd.csv` |
| `ERA5/` | ERA5 reanalysis wind/climate data (`.grib` files) |
| `MEIv2/` | Multivariate ENSO Index v2: `meiv2.data`, `ENSO.csv` |
| `AMO/` | Atlantic Multidecadal Oscillation index: `amo_monthly.txt` |

**Data pipeline scripts** (run in order):

| Script | Description |
|---|---|
| `readERA5.py` | Reads and extracts ERA5 GRIB data *(Python)* |
| `readInterpCTD.R` | Reads and interpolates CTD profiles |
| `readInterpNiskin.R` | Reads and interpolates Niskin bottle data |
| `readInterpPhytoplankton.R` | Reads and interpolates phytoplankton counts |
| `calculateNiskinDepthIntervals.R` | Computes Niskin depth-integrated intervals |
| `calculatePhytoDepthIntervals.R` | Computes phytoplankton depth-integrated intervals |
| `calculatePhytoDiversityFGAbund.R` | Calculates diversity indices and functional group abundances |
| `interpolateData.R` | Interpolates environmental and biological variables |
| `mergeAllData.R` | Merges all processed data into final combined dataset |

**Processed outputs** (`data/processed/`): Cleaned `.RDS` and `.csv` files used directly by plotting scripts, including interpolated phytoplankton counts, Niskin and CTD depth intervals, environmental combined datasets, and genus-level taxonomic information.

---

### `plots/`

R scripts for all main figures, supplemental figures, and review response analyses.

**Main figure scripts:**

| Script | Description |
|---|---|
| `Figure1_Map.R` | Generates study area map |
| `Figure2_ZScoreTimeSeries.R` | All variables z-score time series plot |
| `Figure3_DiversityChlorophyllDepthTimeSeries.R` | Diversity, chlorophyll, and depth time series |
| `Figure4_ClusteringNMDS.R` | Full clustering and NMDS analysis |
| `Figure5_GradientForest.R` | Final gradient forest analysis including time lags |
| `Figure5_GradientForest_preliminaryRuns.R` | Preliminary gradient forest runs |
| `FigureA1_CorrelationClustPlot.R` | Correlation plot (Appendix) |

**`Figure4_Subplots/`** — Individual subfigure scripts for Figure 4:

| Script | Description |
|---|---|
| `Figure4_a_b_YearlyClusterTurnover.R` | Yearly clustering and community turnover |
| `Figure4_a_b_YearlyClusterTurnover_revResponse.R` | Revised version for review response |
| `Figure4_c_NMDS.R` | NMDS ordination plot |
| `Figure4_c_NMDS_EarlyVSLateC1VSC2.R` | NMDS comparing early vs. late cluster periods |
| `Figure4_c_NMDS_seasoncolor.R` | NMDS colored by season |
| `Figure4_d_DensityDistClusters.R` | Cluster-based density distributions |

**`Figure5_GradientForest_*` scripts** — Extended gradient forest analyses for review response, including seasonal analyses, secondary upwelling runs, and supplemental outputs.

**`ReviewResponse/`** — Scripts and outputs generated in response to peer review:

| Script | Description |
|---|---|
| `CommunityClustering_statisticalchecks.R` | Statistical validation of clustering |
| `ClusterComp_WilcoxTable_FINAL.R` | Wilcoxon test comparisons between clusters |
| `EarlyVSLateC1VSC2_statisticstable.R` | Early vs. late period cluster statistics |
| `CalculateProductivity.R` | Primary productivity calculations |
| `ChlorophyllVSPrimaryProductivity.R` | Chlorophyll vs. productivity comparison |
| `Cyano_statisticalchecks.R` | Cyanobacteria-specific checks |
| `PP_statisticalchecks.R` | Primary productivity statistical checks |
| `MissingData_overview.R` / `Final_MissingData_Table_Supplement.R` | Missing data summaries |

**`exports/`** — Generated PDF and `.ai` figures, versioned across manuscript revisions and suitable for further processing in Adobe Illustrator.

---

## Requirements

### R

This project was developed using **R (version 4.3.3)**. Required packages:

```r
install.packages(c(
  "codyn", "corrplot", "cowplot", "data.table", "dendextend", "egg",
  "ggpubr", "ggrepel", "ggspatial", "gradientForest", "kableExtra",
  "mapsf", "marmap", "oce", "patchwork", "RColorBrewer",
  "rnaturalearth", "sf", "tidyverse", "vegan", "viridis", "worrms"
))
```

> Note: Some packages may have additional system dependencies. `gradientForest` may require installation from source or a specific repository.

### Python

One script in the data pipeline uses Python:

- `data/readERA5.py` — requires the `cfgrib` / `eccodes` library for reading ERA5 GRIB files.

```bash
pip install cfgrib
```

---

## Data Sources

The raw phytoplankton, CTD, and Niskin data originate from the **CARIACO Ocean Time-Series Program**, available via [BCO-DMO](https://www.bco-dmo.org/project/2047). Climate indices (ENSO MEIv2, AMO) and ERA5 reanalysis data are publicly available from NOAA and the Copernicus Climate Data Store, respectively.

---

## Citation

If you use or reference this code, please cite the associated Zenodo release:

> Post, B. (2026). *PhytoCariaco* (v1.2). Zenodo. https://doi.org/10.5281/zenodo.18942210

---

## License

This repository is licensed under the [MIT License](LICENSE).

---

## Questions or Contributions

Feedback, bug reports, or contributions are very welcome! Please open an issue or submit a pull request.
