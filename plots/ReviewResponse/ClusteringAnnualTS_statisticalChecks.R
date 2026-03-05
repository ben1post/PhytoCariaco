# ============================================
# ANNUAL COMMUNITY CLUSTERING
# Binary (Jaccard) ONLY - Presence/Absence
# ============================================

library(tidyverse)
library(viridis)
library(cowplot)
library(ggdendro)
library(dendextend)
library(vegan)

# --- SWITCH: Core months ---
USE_CORE_MONTHS <- FALSE  # Set to TRUE to restrict to core months
CORE_MONTHS <- c(1, 2, 4, 5, 9, 11, 12)

# --- 1. LOAD DATA ---
# Ensure the path is correct for your local machine
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# --- 2. PREP DATA ---
ds_genus_yearly <- phyto_counts %>%
  mutate(
    month = as.numeric(format(date, "%m")),
    year = format(date, "%Y")
  ) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species")

# --- 3. APPLY CORE MONTHS FILTER ---
if (USE_CORE_MONTHS) {
  ds_genus_yearly <- ds_genus_yearly %>% filter(month %in% CORE_MONTHS)
  cat("=== USING CORE MONTHS ONLY ===\n")
  cat("Months included:", paste(month.abb[CORE_MONTHS], collapse = ", "), "\n\n")
} else {
  cat("=== USING ALL AVAILABLE MONTHS ===\n\n")
}

# --- 4. CHECK DATA COVERAGE BY YEAR ---
# This is crucial for Jaccard: years with very few months will look "dissimilar"
# simply because we missed species that were actually there.
coverage <- ds_genus_yearly %>%
  group_by(year) %>%
  summarise(
    n_months = n_distinct(month),
    n_genera = n_distinct(Genus),
    .groups = "drop"
  )

cat("=== DATA COVERAGE BY YEAR (Check for low n_months) ===\n")
print(coverage, n = 25)
cat("\n")

# Recommendation: Consider filtering out years with < 5 months if they exist
# ds_genus_yearly <- ds_genus_yearly %>% filter(year %in% coverage$year[coverage$n_months > 4])

# --- 5. AGGREGATE BY GENUS AND YEAR ---
# Using sum() to match original script (Jaccard with binary=TRUE converts to presence/absence anyway)
ds_genus_yearly <- ds_genus_yearly %>%
  group_by(Genus, year) %>%
  summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  arrange(year)

# Filter years (same as original: 1996-2016)
ds_genus_yearly <- ds_genus_yearly %>% filter(year >= 1996 & year <= 2016)

# --- 6. CREATE WIDE MATRIX ---
Mesh_genus_yearly <- pivot_wider(
  ds_genus_yearly,
  names_from = Genus,
  values_from = Total,
  values_fn = function(x) sum(x, na.rm = TRUE),
  values_fill = 0.0
)

cat("=== COMMUNITY MATRIX DIMENSIONS ===\n")
cat(sprintf("Years: %d, Genera: %d\n\n", nrow(Mesh_genus_yearly), ncol(Mesh_genus_yearly) - 1))

# --- 7. CALCULATE JACCARD DISTANCE ---
# binary = TRUE converts abundance to 1/0 automatically
Jaccard_dist <- vegdist(Mesh_genus_yearly[-1], method = "jaccard", binary = TRUE)

Jaccard_dist <- vegdist(Mesh_genus_yearly[-1], method = "bray", binary = FALSE)
Jaccard_mat <- as.matrix(Jaccard_dist)
colnames(Jaccard_mat) <- Mesh_genus_yearly$year
rownames(Jaccard_mat) <- Mesh_genus_yearly$year

# --- 8. HIERARCHICAL CLUSTERING ---
hc_jaccard <- hclust(as.dist(Jaccard_mat)) # default method (complete linkage)
dend_jaccard <- as.dendrogram(hc_jaccard)

# Cut tree (k=2) for coloring
k_jaccard <- cutree(dend_jaccard, k = 2, order_clusters_as_data = FALSE)

# --- 9. PLOTTING SETUP ---
# Color mapping (matches original: cluster 1 = Blue, cluster 2 = Red)
nm_jaccard <- setNames(c("Blue", "Red"), c("1", "2"))
cols_jaccard <- nm_jaccard[as.character(k_jaccard[Mesh_genus_yearly$year])]

# Export cluster year data for use in coloring other plots:
# saveRDS(k_jaccard, "plots/Figure4_Subplots/k_3_coremonths.RDS")

# Prepare Dendrogram object
DENDplot_jaccard <- dend_jaccard %>%
  set("labels_cex", 1.0) %>%
  hang.dendrogram(hang_height = 0.1) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", cols_jaccard, order_value = TRUE) %>%
  set("leaves_bg", cols_jaccard, order_value = TRUE) %>%
  set("branches_lwd", 0.5)

# --- 10. GENERATE PLOT ---
# Matching original plot settings
DendroPLOT <- ggplot(DENDplot_jaccard, horiz = FALSE, offset_labels = -0.03) +
  theme_cowplot(font_size = 20) +
  ylab("Height") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.clip = "off"
  ) +
  scale_y_continuous(limits = c(0., 1.0)) +
  scale_x_reverse()

print(DendroPLOT)

# --- 11. EXPORT ---
# ggsave("plots/Jaccard_Annual_Clustering_CoreMonths.png", DendroPLOT, width = 10, height = 6, dpi = 300)