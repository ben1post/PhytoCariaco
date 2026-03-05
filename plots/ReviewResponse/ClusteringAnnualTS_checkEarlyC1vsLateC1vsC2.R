# ============================================
# PERIOD-LEVEL COMMUNITY CLUSTERING
# Aggregate by Early C1, C2, Late C1
# Binary (Jaccard) and Dominance (Bray-Curtis)
# ============================================

library(tidyverse)
library(viridis)
library(cowplot)
library(ggdendro)
library(dendextend)
library(ggpubr)
library(egg)
library(vegan)

# --- SWITCH: Core months or all months ---
USE_CORE_MONTHS <- FALSE  # Set to TRUE to restrict to core months
CORE_MONTHS <- c(1, 2, 4, 5, 9, 11, 12)

# --- PERIOD DEFINITIONS ---
# Based on manuscript: Early C1 (1996-2004), C2 (2005-2013), Late C1 (2014-2016)
PERIOD_BREAKS <- list(
  Early_C1 = 1996:2004,
  C2 = 2005:2013,
  Late_C1 = 2014:2016
)

# --- 1. LOAD DATA ---
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# --- 2. PREP DATA ---
ds_genus <- phyto_counts %>%
  mutate(
    month = as.numeric(format(date, "%m")),
    year = as.numeric(format(date, "%Y"))
  ) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species")

# --- 3. APPLY CORE MONTHS FILTER (if enabled) ---
if (USE_CORE_MONTHS) {
  ds_genus <- ds_genus %>% filter(month %in% CORE_MONTHS)
  cat("=== USING CORE MONTHS ONLY ===\n")
  cat("Months included:", paste(month.abb[CORE_MONTHS], collapse = ", "), "\n\n")
} else {
  cat("=== USING ALL AVAILABLE MONTHS ===\n\n")
}

# --- 4. ASSIGN PERIODS ---
ds_genus <- ds_genus %>%
  mutate(
    period = case_when(
      year %in% PERIOD_BREAKS$Early_C1 ~ "Early_C1",
      year %in% PERIOD_BREAKS$C2 ~ "C2",
      year %in% PERIOD_BREAKS$Late_C1 ~ "Late_C1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period))

# --- 5. CHECK DATA COVERAGE BY PERIOD ---
coverage <- ds_genus %>%
  group_by(period) %>%
  summarise(
    year_range = paste(min(year), max(year), sep = "-"),
    n_months = n_distinct(paste(year, month)),
    n_genera = n_distinct(Genus),
    .groups = "drop"
  ) %>%
  mutate(period = factor(period, levels = c("Early_C1", "C2", "Late_C1"))) %>%
  arrange(period)

cat("=== DATA COVERAGE BY PERIOD ===\n")
print(coverage)
cat("\n")

# # --- 6. AGGREGATE BY GENUS AND PERIOD ---
# ds_genus_period <- ds_genus %>%
#   group_by(Genus, period) %>%
#   summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop")



# --- 6. AGGREGATE BY GENUS AND PERIOD (normalized by sampling effort) ---

# First, get number of monthly samples per period
n_months_per_period <- ds_genus %>%
  group_by(period) %>%
  summarise(n_samples = n_distinct(paste(year, month)), .groups = "drop")

ds_genus_period <- ds_genus %>%
  group_by(Genus, period) %>%
  summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  left_join(n_months_per_period, by = "period") %>%
  mutate(Total = Total / n_samples) %>%
  select(-n_samples)



# --- 7. CREATE WIDE MATRIX ---
Mesh_genus_period <- pivot_wider(
  ds_genus_period,
  names_from = Genus,
  values_from = Total,
  values_fn = function(x) sum(x, na.rm = TRUE),
  values_fill = 0.0
)

# Ensure correct row order
Mesh_genus_period <- Mesh_genus_period %>%
  mutate(period = factor(period, levels = c("Early_C1", "C2", "Late_C1"))) %>%
  arrange(period)

cat("=== COMMUNITY MATRIX DIMENSIONS ===\n")
cat(sprintf("Periods: %d, Genera: %d\n\n", nrow(Mesh_genus_period), ncol(Mesh_genus_period) - 1))

# --- 8. CALCULATE DISTANCE MATRICES ---

# Binary (Jaccard) - presence/absence
Jaccard_dist <- vegdist(Mesh_genus_period[-1], method = "jaccard", binary = TRUE)
Jaccard_mat <- as.matrix(Jaccard_dist)
colnames(Jaccard_mat) <- Mesh_genus_period$period
rownames(Jaccard_mat) <- Mesh_genus_period$period

# Dominance (Bray-Curtis) - abundance-weighted
BrayCurtis_dist <- vegdist(Mesh_genus_period[-1], method = "bray", binary = FALSE)
BrayCurtis_mat <- as.matrix(BrayCurtis_dist)
colnames(BrayCurtis_mat) <- Mesh_genus_period$period
rownames(BrayCurtis_mat) <- Mesh_genus_period$period

# --- 9. REPORT DISTANCE/SIMILARITY MATRICES ---
cat("=== JACCARD DISTANCE (Binary) ===\n")
print(round(Jaccard_mat, 3))
cat("\n=== JACCARD SIMILARITY (1 - distance) ===\n")
print(round(1 - Jaccard_mat, 3))

cat("\n=== BRAY-CURTIS DISTANCE (Abundance) ===\n")
print(round(BrayCurtis_mat, 3))
cat("\n=== BRAY-CURTIS SIMILARITY (1 - distance) ===\n")
print(round(1 - BrayCurtis_mat, 3))

# --- 10. HIERARCHICAL CLUSTERING ---

# Jaccard clustering
hc_jaccard <- hclust(as.dist(Jaccard_mat), method = "ward.D2")
dend_jaccard <- as.dendrogram(hc_jaccard)

# Bray-Curtis clustering
hc_bray <- hclust(as.dist(BrayCurtis_mat), method = "ward.D2")
dend_bray <- as.dendrogram(hc_bray)

# --- 11. REPORT CLUSTERING ORDER ---
cat("\n=== CLUSTERING RESULTS ===\n")
cat("Jaccard clustering order:", paste(hc_jaccard$labels[hc_jaccard$order], collapse = " <-> "), "\n")
cat("Bray-Curtis clustering order:", paste(hc_bray$labels[hc_bray$order], collapse = " <-> "), "\n")

# Identify which periods cluster together first (lowest merge height)
cat("\n=== FIRST MERGE (most similar pair) ===\n")
cat("Jaccard: ")
if (hc_jaccard$merge[1, 1] < 0 && hc_jaccard$merge[1, 2] < 0) {
  cat(hc_jaccard$labels[-hc_jaccard$merge[1, 1]], "and", 
      hc_jaccard$labels[-hc_jaccard$merge[1, 2]], "\n")
}
cat("Bray-Curtis: ")
if (hc_bray$merge[1, 1] < 0 && hc_bray$merge[1, 2] < 0) {
  cat(hc_bray$labels[-hc_bray$merge[1, 1]], "and", 
      hc_bray$labels[-hc_bray$merge[1, 2]], "\n")
}

# --- 12. PREPARE DENDROGRAMS ---

# Color by period
period_colors <- c("Early_C1" = "#1f77b4", "C2" = "#d62728", "Late_C1" = "#2ca02c")

# Jaccard dendrogram
DENDplot_jaccard <- dend_jaccard %>%
  set("labels_cex", 1.2) %>%
  hang.dendrogram(hang_height = 0.02) %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 2) %>%
  set("branches_lwd", 1)

# Bray-Curtis dendrogram
DENDplot_bray <- dend_bray %>%
  set("labels_cex", 1.2) %>%
  hang.dendrogram(hang_height = 0.02) %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 2) %>%
  set("branches_lwd", 1)

# --- 13. PLOT DENDROGRAMS ---

# Jaccard plot
DendroPLOT_jaccard <- ggplot(DENDplot_jaccard, horiz = FALSE, offset_labels = -0.01) +
  theme_cowplot(font_size = 16) +
  ylab("Height") +
  ggtitle("Jaccard (Presence/Absence)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  scale_x_reverse()

# Bray-Curtis plot
DendroPLOT_bray <- ggplot(DENDplot_bray, horiz = FALSE, offset_labels = -0.01) +
  theme_cowplot(font_size = 16) +
  ylab("Height") +
  ggtitle("Bray-Curtis (Abundance)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  scale_x_reverse()

# Combined plot
combined_plot <- plot_grid(
  DendroPLOT_jaccard, DendroPLOT_bray,
  ncol = 2, labels = c("A", "B")
)

print(combined_plot)

w# --- 14. INTERPRETATION SUMMARY ---
cat("\n=== INTERPRETATION ===\n")
cat("If Early_C1 and Late_C1 cluster together first: supports 'return' narrative\n")
cat("If C2 and Late_C1 cluster together first: suggests community did not return\n")
