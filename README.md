# PhytoCariaco
Analysis of Phytoplankton Community Changes and Bottom Up Drivers in the Cariaco Basin Venezuela

# Phytoplankton Community Dynamics in the Cariaco Basin

This repository contains the R code used for the data analysis and figure generation presented in the manuscript:

**"Decadal-Scale Oscillations in Phytoplankton Community Composition in the Cariaco Basin, Venezuela"**  
*Submitted to [Journal Name], 2025*

## Overview

The analysis investigates long-term changes in phytoplankton community structure in the Cariaco Basin using presence-absence and abundance data, hierarchical clustering, and ordination methods. The code in this repository reproduces all key results, including:

- Community clustering based on Jaccard and Bray-Curtis distances  
- Non-metric multidimensional scaling (NMDS) of phytoplankton abundance data  
- Environmental vector fitting to community ordination  
- Calculation of diversity and turnover metrics  
- All figures included in the manuscript

## Repository Structure

📁 data/ # Cleaned data used in analysis (not raw data)
📁 scripts/ # R scripts for analysis and plotting
├── 01_preprocessing.R
├── 02_clustering.R
├── 03_nmds_analysis.R
├── 04_diversity_metrics.R
├── 05_figure_generation.R
📁 output/ # Generated figures and intermediate results
📄 README.md # This file
📄 CITATION.cff # Citation information for this repository


## Requirements

This project was developed using R (version 4.3.2). Required packages include:

- `vegan`
- `tidyverse`
- `cluster`
- `ggplot2`
- `patchwork`
- `RColorBrewer`

You can install them using:

```R
install.packages(c("vegan", "tidyverse", "cluster", "ggplot2", "patchwork", "RColorBrewer"))
