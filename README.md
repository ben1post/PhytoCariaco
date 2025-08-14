# PhytoCariaco: Phytoplankton Community Dynamics in the Cariaco Basin

This repository contains the R code used for the data analysis and figure generation presented in the manuscript:

**"Reversible Regime Change: Climate-driven phytoplankton community shifts in the Cariaco Basin, Venezuela"**  
*Submitted to JGR Biogeosciences, 2025*
Authored by Benjamin Post, Acevedo-Trejos, Chakraborty, Barton, and Merico.

All code was written by Benjamin Post.


## Overview

The analysis investigates long-term changes in phytoplankton community structure in the Cariaco Basin using time series analysis and statistics, hierarchical clustering, ordination methods and an analysis of potential drivers of change using the gradient forest algorithm.

## Repository Structure

📁 data/ # Raw data, data pipeline and cleaned data used in analysis
📁 plots/ # R scripts for analysis and plotting
├── Figure1_Map.R
├── Figure2_ZScoreTimeSeries.R
├── Figure3_DiversityChlorophyllDepthTimeSeries.R
├── Figure4_ClusteringNMDS.R
├── 📁 Figure4_Subplots/ # Scripts to generate subfigures of Figure 4
├──├── Figure4_a_b_YearlyClusterTurnover.R
├──├── Figure4_c_NMDS.R
├──├── Figure4_d_DensityDistClusters
├── Figure5_GradientForest_preliminaryRuns.R
├── FigureA1_CorrelationClustPlot.R
├── 📁 export/ # Generated figures which were further processed with Adobe Illustrator
📄 README.md # This file
📄 LICENSE # License information for this repository


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
Note: Some of these packages might have additional dependencies, based on you system.

## Questions or Contributions

Feedback, bug reports, or contributions are very welcome! Please open an issue or submit a pull request.