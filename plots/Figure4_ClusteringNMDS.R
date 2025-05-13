library(tidyverse)

require(cowplot)
library(ggdendro)
library(dendextend)

library(cowplot)
library(ggpubr)
library(egg)

library(codyn)


# read interpolated phytoplankton cell count data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

unique(phyto_counts$FuncGroup)


source("plots/Figure4_Subplots/Figure4_a_b_YearlyClusterTurnover.R")

source("plots/Figure4_Subplots/Figure4_c_NMDS.R")


###### COMBINE INTO SINGLE PLOT #######


options(repr.plot.width=13, repr.plot.height=7)

DendTurnPLOT <- plot_grid(DendroPLOT, TurnoverPLOT, ncol=1, align="hv", labels=c("a","b"))
DendTurnPLOT

CLUSTERPLOT <- plot_grid(DendTurnPLOT,NMDSplot, align="v", rel_widths=c(1,1), labels = c('','c'))#, label_size = 22)
CLUSTERPLOT


source("plots/Figure4_Subplots/Figure4_d_DensityDistClusters.R")


###### EXPORT PLOT!
library("patchwork")

options(repr.plot.width=12, repr.plot.height=14)

FullNewPlotGrid <- DendroPLOT + TurnoverPLOT + NMDSplot + legenddd +
  plot_layout(
    design="
    AACC
    AACC
    BBCC
    BBCC
    DDDD
    ", heights=c(1,1,1,1,0.5)) + plot_annotation(tag_levels="a") +
  theme(plot.tag = element_text(size = 20))

FullNewPlotGrid

ggsave("plots/exports/Figure4_ClusteringNMDS_v2.pdf",FullNewPlotGrid, width=12, height=7)


BottomPlotGrid <- A + B + C + D + E + F + G + H + I + J + K + L + M + N + O +
  plot_layout(
  design="
    ABCDE
    FGHIJ
    KLMNO") + plot_annotation(tag_levels="a") +
  theme(plot.tag = element_text(size = 20))
  
BottomPlotGrid

ggsave("plots/exports/Figure5_ClustCompPlot_v2.pdf",BottomPlotGrid, width=13, height=7)
