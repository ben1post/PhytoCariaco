# PhytoCariaco: Phytoplankton Community Dynamics in the Cariaco Basin

This repository contains the R code used for the data analysis and figure generation presented in the manuscript:

**"Reversible Regime Change: Climate-driven phytoplankton community shifts in the Cariaco Basin, Venezuela"**  
*Submitted to JGR Biogeosciences, 2025*

Authored by Post, Acevedo-Trejos, Chakraborty, Barton, and Merico.

All code was written by Post.


## Overview

The analysis investigates long-term changes in phytoplankton community structure in the Cariaco Basin using time series analysis and statistics, hierarchical clustering, ordination methods and an analysis of potential drivers of change using the gradient forest algorithm.

## Repository Structure

- `data/`  
  Raw data, data pipeline, and cleaned data used in analysis

- `plots/`  
  R scripts for analysis and plotting  
  - `Figure1_Map.R` — Generates study area map  
  - `Figure2_ZScoreTimeSeries.R` — All variables z-score time series plot
  - `Figure3_DiversityChlorophyllDepthTimeSeries.R` — Diversity, chlorophyll, and depth time series  
  - `Figure4_ClusteringNMDS.R` — Clustering and NMDS analysis  
  - `Figure4_Subplots/`  
    Scripts for subfigures of Figure 4:  
    - `Figure4_a_b_YearlyClusterTurnover.R` — Yearly clustering and turnover  
    - `Figure4_c_NMDS.R` — NMDS ordination plot  
    - `Figure4_d_DensityDistClusters` — Cluster-based density distributions  
  - `Figure5_GradientForest_preliminaryRuns.R` — Preliminary gradient forest analysis  
  - `Figure5_GradientForest.R` — Final gradient forest analysis incl. time lags 
  - `FigureA1_CorrelationClustPlot.R` — Correlation plot in Appendix

- `export/`  
  Generated figures which were further processed with Adobe Illustrator

- `README.md`  
  This file

- `LICENSE`  
  License information for this repository



## Requirements

This project was developed using R (version 4.3.3). Required packages include:

- codyn  
- corrplot  
- cowplot  
- data.table  
- dendextend  
- egg  
- ggplot2 (included in tidyverse)  
- ggpubr  
- ggrepel  
- ggspatial  
- gradientForest  
- kableExtra  
- mapsf  
- marmap  
- oce  
- patchwork  
- RColorBrewer  
- rnaturalearth  
- sf  
- tidyverse  
- vegan  
- viridis  
- worrms  

You can install them using:

```R
install.packages(c(
  "codyn", "corrplot", "cowplot", "data.table", "dendextend", "egg",
  "ggpubr", "ggrepel", "ggspatial", "gradientForest", "kableExtra",
  "mapsf", "marmap", "oce", "patchwork", "RColorBrewer",
  "rnaturalearth", "sf", "tidyverse", "vegan", "viridis", "worrms"
))
```
Note: Some of these packages might have additional dependencies, based on your system.

## Questions or Contributions

Feedback, bug reports, or contributions are very welcome! Please open an issue or submit a pull request.