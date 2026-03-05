# =============================================================================
# PHYTOPLANKTON COMMUNITY CLUSTERING ANALYSIS
# Jaccard (presence/absence) and Bray-Curtis (abundance) comparison
# =============================================================================

library(tidyverse)
library(vegan)
library(dendextend)
library(cowplot)

phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# =============================================================================
# PART 1: JACCARD CLUSTERING (Presence/Absence)
# =============================================================================

# Prepare data - sum counts per genus per year
ds_jac <- phyto_counts %>%
  mutate(year = format(date, "%Y")) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, year) %>%
  summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  filter(year >= 1996 & year <= 2016)

mat_jac <- pivot_wider(ds_jac, names_from = Genus, values_from = Total,
                       values_fn = function(x) sum(x, na.rm = TRUE), 
                       values_fill = 0.0)

# Jaccard clustering
jac <- vegdist(mat_jac[-1], method = "jaccard", binary = TRUE)
jac_mat <- as.matrix(jac)
colnames(jac_mat) <- mat_jac$year
rownames(jac_mat) <- mat_jac$year

hc_jac <- hclust(as.dist(jac_mat))
dend_jac <- as.dendrogram(hc_jac)
k_jac <- cutree(dend_jac, k = 3, order_clusters_as_data = FALSE)

# Colors for Jaccard clusters
cluster_cols <- c("1" = "#B2182B", "2" = "#92C5DE", "3" = "#2166AC")
cols_jac <- cluster_cols[as.character(k_jac[mat_jac$year])]

# Plot Jaccard dendrogram
dend_jac_styled <- dend_jac %>%
  set("labels_cex", 1.) %>%
  hang.dendrogram(hang_height = 0.1) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", cols_jac, order_value = TRUE) %>%
  set("leaves_bg", cols_jac, order_value = TRUE) %>%
  set("branches_lwd", 0.5)

plot_jac <- ggplot(dend_jac_styled, horiz = FALSE, offset_labels = -0.03) +
  theme_cowplot(font_size = 20) +
  ylab("Height (Jaccard)") +
  theme(
    axis.title.x = element_blank(),
    strip.clip = "off"
  ) +
  scale_y_continuous(limits = c(0., 0.7)) +
  scale_x_reverse()

print(plot_jac)

# =============================================================================
# PART 2: BRAY-CURTIS CLUSTERING (Abundance-based)
# =============================================================================

# Prepare data - use MEAN to account for variable months sampled per year
ds_bc <- phyto_counts %>%
  mutate(year = format(date, "%Y")) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, year) %>%
  summarise(MeanCount = mean(counts, na.rm = TRUE), .groups = "drop") %>%
  filter(year >= 1996 & year <= 2016)

mat_bc <- pivot_wider(ds_bc, names_from = Genus, values_from = MeanCount,
                      values_fn = function(x) mean(x, na.rm = TRUE), 
                      values_fill = 0.0)

# Hellinger transformation
mat_bc_transformed <- decostand(mat_bc[-1], method = "hellinger")

# Bray-Curtis clustering
bc <- vegdist(mat_bc_transformed, method = "bray")
bc_mat <- as.matrix(bc)
colnames(bc_mat) <- mat_bc$year
rownames(bc_mat) <- mat_bc$year

hc_bc <- hclust(as.dist(bc_mat))
dend_bc <- as.dendrogram(hc_bc)
k_bc <- cutree(dend_bc, k = 3, order_clusters_as_data = FALSE)

# Use Jaccard colors for comparison
cols_bc <- cols_jac

# Plot Bray-Curtis dendrogram
dend_bc_styled <- dend_bc %>%
  set("labels_cex", 1.) %>%
  hang.dendrogram(hang_height = 0.1) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", cols_bc, order_value = TRUE) %>%
  set("leaves_bg", cols_bc, order_value = TRUE) %>%
  set("branches_lwd", 0.5)

plot_bc <- ggplot(dend_bc_styled, horiz = FALSE, offset_labels = -0.03) +
  theme_cowplot(font_size = 20) +
  ylab("Height (Bray-Curtis)") +
  theme(
    axis.title.x = element_blank(),
    strip.clip = "off"
  ) +
  scale_y_continuous(limits = c(0., 0.8)) +
  scale_x_reverse()

print(plot_bc)

# =============================================================================
# PART 3: COMPARE DISTANCE MATRICES
# =============================================================================

# Mantel test
mantel_result <- mantel(jac, bc, method = "spearman")
print(mantel_result)

# Heatmaps
par(mfrow = c(1, 2))
heatmap(as.matrix(jac), main = "Jaccard", symm = TRUE)
heatmap(as.matrix(bc), main = "Bray-Curtis", symm = TRUE)
par(mfrow = c(1, 1))

# =============================================================================
# PART 4: RICHNESS VS ABUNDANCE COMPARISON
# =============================================================================

yearly_summary <- phyto_counts %>%
  mutate(year = as.numeric(format(date, "%Y"))) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  filter(year >= 1996 & year <= 2016) %>%
  filter(counts > 0) %>%
  group_by(year) %>%
  summarise(
    richness = n_distinct(Genus),
    total_abundance = sum(counts, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    richness_scaled = (richness - mean(richness)) / sd(richness),
    abundance_scaled = (total_abundance - mean(total_abundance)) / sd(total_abundance)
  ) %>%
  pivot_longer(cols = c(richness_scaled, abundance_scaled), 
               names_to = "metric", values_to = "value")

plot_richness_abundance <- ggplot(yearly_summary, aes(x = year, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("richness_scaled" = "darkgreen", "abundance_scaled" = "purple"),
                     labels = c("Abundance (scaled)", "Richness (scaled)")) +
  theme_cowplot() +
  labs(y = "Z-score", x = "Year", color = "")

print(plot_richness_abundance)

# =============================================================================
# PART 5: TAXA SUMMARY
# =============================================================================

taxa_summary <- ds_jac %>%
  group_by(Genus) %>%
  summarise(
    occurrence = n_distinct(year),
    mean_abundance = mean(Total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    type = case_when(
      occurrence > 15 & mean_abundance > median(mean_abundance) ~ "Core dominant",
      occurrence > 15 ~ "Core rare",
      mean_abundance > median(mean_abundance) ~ "Sporadic dominant",
      TRUE ~ "Sporadic rare"
    )
  )

print(table(taxa_summary$type))

# =============================================================================
# PART 6: PERMANOVA
# =============================================================================

# --- PERMANOVA for Jaccard (3 clusters) ---
cluster_groups_jac <- factor(k_jac[mat_jac$year])

set.seed(42)
permanova_jac <- adonis2(
  mat_jac[-1] ~ cluster_groups_jac,
  method = "jaccard",
  binary = TRUE,
  permutations = 999
)
cat("\n=== PERMANOVA: Jaccard (3 clusters) ===\n")
print(permanova_jac)

# Test homogeneity of dispersion
betadisper_jac <- betadisper(jac, cluster_groups_jac)
betadisper_jac_test <- permutest(betadisper_jac, permutations = 999)
cat("\n=== Dispersion test: Jaccard ===\n")
print(betadisper_jac_test)

# --- PERMANOVA for Bray-Curtis (3 clusters) ---
cluster_groups_bc <- factor(k_bc[mat_bc$year])

set.seed(42)
permanova_bc <- adonis2(
  mat_bc_transformed ~ cluster_groups_bc,
  method = "bray",
  permutations = 999
)
cat("\n=== PERMANOVA: Bray-Curtis (3 clusters) ===\n")
print(permanova_bc)

# Test homogeneity of dispersion
betadisper_bc <- betadisper(bc, cluster_groups_bc)
betadisper_bc_test <- permutest(betadisper_bc, permutations = 999)
cat("\n=== Dispersion test: Bray-Curtis ===\n")
print(betadisper_bc_test)

# =============================================================================
# PART 7: SIMPER ANALYSIS
# =============================================================================

# Run SIMPER for Bray-Curtis clusters
simper_bc <- simper(mat_bc_transformed, group = k_bc, permutations = 999)

cat("\n=== SIMPER Summary ===\n")
summary(simper_bc)

# --- Function to summarize SIMPER results (CORRECTED) ---
summarize_simper <- function(simper_obj, p_threshold = 0.05) {
  
  comparisons <- names(simper_obj)
  
  results_list <- lapply(comparisons, function(comp) {
    sim <- simper_obj[[comp]]
    
    # Use names() instead of rownames() - sim$average is a named vector
    df <- data.frame(
      Genus = names(sim$average),
      Comparison = comp,
      Contribution = as.numeric(sim$average),
      SD = as.numeric(sim$sd),
      Ratio = as.numeric(sim$average / sim$sd),
      Abund_GroupA = as.numeric(sim$ava),
      Abund_GroupB = as.numeric(sim$avb),
      Cumsum = as.numeric(sim$cusum),
      p_value = as.numeric(sim$p),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
    # Determine which group the taxon indicates
    groups <- strsplit(comp, "_")[[1]]
    df$Indicates <- ifelse(df$Abund_GroupA > df$Abund_GroupB, groups[1], groups[2])
    
    # Fold difference
    df$FoldDiff <- pmax(df$Abund_GroupA, df$Abund_GroupB) / 
      pmax(pmin(df$Abund_GroupA, df$Abund_GroupB), 0.001)
    
    return(df)
  })
  
  results_df <- do.call(rbind, results_list)
  
  # Significant indicators
  sig_indicators <- results_df %>%
    filter(p_value <= p_threshold) %>%
    arrange(Comparison, desc(Contribution))
  
  return(list(
    all_results = results_df,
    significant = sig_indicators
  ))
}

# Re-run with corrected function
simper_summary <- summarize_simper(simper_bc, p_threshold = 0.05)

cat("\n=== Significant Indicator Taxa (p <= 0.05) ===\n")
print(simper_summary$significant %>% 
        select(Comparison, Genus, Indicates, Contribution, Ratio, 
               Abund_GroupA, Abund_GroupB, p_value) %>%
        group_by(Comparison) %>%
        slice_head(n = 10),
      n = 10)

# --- Summary table: Top indicators per cluster ---
create_cluster_indicator_table <- function(simper_summary, k_clusters) {
  
  cluster_names <- sort(unique(k_clusters))
  
  summary_list <- lapply(cluster_names, function(clust) {
    
    indicators <- simper_summary$significant %>%
      filter(Indicates == as.character(clust)) %>%
      group_by(Genus) %>%
      summarise(
        mean_contribution = mean(Contribution),
        mean_abundance = mean(pmax(Abund_GroupA, Abund_GroupB)),
        n_comparisons_significant = n(),
        mean_ratio = mean(Ratio, na.rm = TRUE),
        min_p = min(p_value),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_contribution)) %>%
      slice_head(n = 10)
    
    indicators$Cluster <- clust
    return(indicators)
  })
  
  do.call(rbind, summary_list)
}

cluster_indicators <- create_cluster_indicator_table(simper_summary, k_bc)

cat("\n=== Top Indicator Taxa per Cluster ===\n")
print(cluster_indicators %>% 
        select(Cluster, Genus, mean_contribution, mean_abundance, 
               n_comparisons_significant, mean_ratio, min_p) %>%
        arrange(Cluster, desc(mean_contribution)),
      n = 30)



