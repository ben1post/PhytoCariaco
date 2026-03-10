# =============================================================================
# PHYTOPLANKTON COMMUNITY CLUSTERING ANALYSIS
# Comparison of Jaccard vs Bray-Curtis, Complete vs Ward's method
# =============================================================================
library(tidyverse)
library(vegan)
library(dendextend)
library(cowplot)

phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# =============================================================================
# PREPARE DATA
# =============================================================================

# Jaccard data (presence/absence) - use sum
ds_jac <- phyto_counts %>%
  mutate(year = format(date, "%Y")) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, year) %>%
  summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  filter(year >= 1996 & year <= 2016)

mat_jac <- pivot_wider(ds_jac, names_from = Genus, values_from = Total,
                       values_fn = function(x) sum(x, na.rm = TRUE), 
                       values_fill = 0.0)

# Bray-Curtis data (abundance) - use mean, then cube-root transform
ds_bc <- phyto_counts %>%
  mutate(year = format(date, "%Y")) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, year) %>%
  summarise(MeanCount = mean(counts, na.rm = TRUE), .groups = "drop") %>%
  filter(year >= 1996 & year <= 2016)

mat_bc <- pivot_wider(ds_bc, names_from = Genus, values_from = MeanCount,
                      values_fn = function(x) mean(x, na.rm = TRUE), 
                      values_fill = 0.0)

# Cube-root transformation for Bray-Curtis
mat_bc_transformed <- mat_bc
mat_bc_transformed[-1] <- (mat_bc[-1])^(1/3)

# =============================================================================
# CALCULATE DISTANCE MATRICES
# =============================================================================

# Jaccard distance
jac_dist <- vegdist(mat_jac[-1], method = "jaccard", binary = TRUE)
jac_mat <- as.matrix(jac_dist)
colnames(jac_mat) <- mat_jac$year
rownames(jac_mat) <- mat_jac$year

# Bray-Curtis distance
bc_dist <- vegdist(mat_bc_transformed[-1], method = "bray")
bc_mat <- as.matrix(bc_dist)
colnames(bc_mat) <- mat_bc$year
rownames(bc_mat) <- mat_bc$year

# =============================================================================
# GET REFERENCE CLUSTER ASSIGNMENTS (Jaccard + Complete linkage)
# =============================================================================
hc_reference <- hclust(as.dist(jac_mat), method = "complete")
dend_reference <- as.dendrogram(hc_reference)
k_3_reference <- cutree(dend_reference, k = 3, order_clusters_as_data = FALSE)

# Color mapping:
# cutree 1 = Cluster 2 (dark red)
# cutree 2 = Late Cluster 1 (light blue)
# cutree 3 = Early Cluster 1 (dark blue)
nm2 <- c("1" = "#B2182B", "2" = "#92C5DE", "3" = "#2166AC")

# =============================================================================
# FUNCTION TO CREATE DENDROGRAM WITH CONSISTENT COLORING
# =============================================================================
create_dend_plot <- function(dist_matrix, method, k_3_reference, color_map, y_limit = 0.75) {
  
  
  # Hierarchical clustering with specified method
  hc <- hclust(dist_matrix, method = method)
  dend <- as.dendrogram(hc)
  
  # Get the order of leaves in the dendrogram
  leaf_order <- labels(dend)
  
  # Create color vector in dendrogram leaf order
  cols_ordered <- color_map[as.character(k_3_reference[leaf_order])]
  
  # Style the dendrogram
  dend_styled <- dend %>%
    set("labels_cex", 0.8) %>%
    hang.dendrogram(hang_height = 0.1) %>% 
    set("leaves_pch", 19) %>% 
    set("leaves_col", cols_ordered) %>% 
    set("leaves_cex", 1.5) %>%
    set("branches_lwd", 0.5)
  
  # Create ggplot (no individual titles)
  p <- ggplot(dend_styled, horiz = FALSE, offset_labels = -0.03) + 
    theme_cowplot(font_size = 12) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
    ) + 
    scale_y_continuous(limits = c(0., y_limit)) + 
    scale_x_reverse()
  
  return(p)
}

# =============================================================================
# CREATE ALL FOUR DENDROGRAM PLOTS
# =============================================================================

# A: Jaccard + Complete
plot_jac_complete <- create_dend_plot(
  as.dist(jac_mat), "complete", k_3_reference, nm2, y_limit = 0.75
)

# B: Jaccard + Ward's method
plot_jac_ward <- create_dend_plot(
  as.dist(jac_mat), "ward.D2", k_3_reference, nm2, y_limit = 1.25
)

# C: Bray-Curtis + Complete
plot_bc_complete <- create_dend_plot(
  as.dist(bc_mat), "complete", k_3_reference, nm2, y_limit = 0.85
)

# D: Bray-Curtis + Ward's method
plot_bc_ward <- create_dend_plot(
  as.dist(bc_mat), "ward.D2", k_3_reference, nm2, y_limit = 1.55
)

# =============================================================================
# CREATE COLUMN HEADERS (Linkage methods)
# =============================================================================
header_complete <- ggdraw() + 
  draw_label("Complete linkage", fontface = "bold", size = 12, hjust = 0.5)

header_ward <- ggdraw() + 
  draw_label("Ward's method", fontface = "bold", size = 12, hjust = 0.5)

# =============================================================================
# CREATE ROW LABELS (Distance metrics)
# =============================================================================
label_jaccard <- ggdraw() + 
  draw_label("Jaccard", fontface = "bold", size = 12, angle = 90, hjust = 0.5)

label_braycurtis <- ggdraw() + 
  draw_label("Bray-Curtis", fontface = "bold", size = 12, angle = 90, hjust = 0.5)

# =============================================================================
# ASSEMBLE THE GRID
# =============================================================================

# Create the 2x2 grid of dendrograms
dend_grid <- plot_grid(
  plot_jac_complete, plot_jac_ward,
  plot_bc_complete, plot_bc_ward,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14
)

# Add column headers
header_row <- plot_grid(
  NULL, header_complete, header_ward,
  ncol = 3, rel_widths = c(0.08, 1, 1)
)

# Add row labels to the dendrogram grid
dend_with_row_labels <- plot_grid(
  label_jaccard, 
  plot_grid(plot_jac_complete, plot_jac_ward, ncol = 2, labels = c("A", "B"), label_size = 14),
  label_braycurtis, 
  plot_grid(plot_bc_complete, plot_bc_ward, ncol = 2, labels = c("C", "D"), label_size = 14),
  ncol = 2, nrow = 2,
  rel_widths = c(0.08, 1),
  rel_heights = c(1, 1)
)

# Combine header row with the main grid
combined_plot <- plot_grid(
  header_row,
  dend_with_row_labels,
  ncol = 1,
  rel_heights = c(0.06, 1)
)

# Add overall title
combined_plot_titled <- ggdraw() +
  draw_label("Comparison of distance metrics and linkage methods", 
             x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, size = 14, fontface = "bold") +
  draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 0.94)

combined_plot_titled

# Save
ggsave("plots/exports/Supplemental_ClusteringComparison.pdf", combined_plot_titled, 
       width = 8, height = 6)
