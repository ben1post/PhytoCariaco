# ============================================
# SENSITIVITY ANALYSIS: Testing Clustering Robustness
# Option A: Use months with ≥90% coverage
# Option B: Exclude years with poor sampling
# ============================================

library(tidyverse)
library(vegan)
library(dendextend)
library(cowplot)

# --- 1. LOAD DATA ---
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# --- 2. CREATE MONTHLY PRESENCE MATRIX ---
monthly_coverage <- phyto_counts %>%
  mutate(
    month = as.numeric(format(date, "%m")),
    year = as.numeric(format(date, "%Y"))
  ) %>%
  filter(year >= 1996 & year <= 2016) %>%
  group_by(year, month) %>%
  summarise(sampled = n() > 0, .groups = "drop") %>%
  complete(year = 1996:2016, month = 1:12, fill = list(sampled = FALSE))

# --- 3. DEFINE FILTERING OPTIONS ---

# OPTION A: Months with ≥90% coverage across all years
# From diagnostic: Feb (95.2%), Mar (90.5%), Jun (90.5%), Dec (90.5%)
OPTION_A_MONTHS <- c(2, 3, 6, 12)

# OPTION B: Years with good sampling (≥10 months)
coverage_by_year <- monthly_coverage %>%
  group_by(year) %>%
  summarise(n_months = sum(sampled), .groups = "drop")

OPTION_B_YEARS <- coverage_by_year %>%
  filter(n_months >= 10) %>%
  pull(year)

cat("=== OPTION A: High-coverage months ===\n")
cat("Months:", paste(month.abb[OPTION_A_MONTHS], collapse = ", "), "\n\n")

cat("=== OPTION B: Well-sampled years ===\n")
cat("Years included:", paste(OPTION_B_YEARS, collapse = ", "), "\n")
cat("Years excluded:", paste(setdiff(1995:2017, OPTION_B_YEARS), collapse = ", "), "\n\n")

# --- 4. FUNCTION TO RUN CLUSTERING ---
run_clustering <- function(phyto_data, label) {
  
  # Aggregate by genus and year
  ds <- phyto_data %>%
    mutate(year = format(date, "%Y")) %>%
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
    group_by(Genus, year) %>%
    summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop") %>%
    filter(year >= "1996" & year <= "2016")
  
  # Wide matrix
  mat <- pivot_wider(ds, names_from = Genus, values_from = Total,
                     values_fn = function(x) sum(x, na.rm = TRUE), 
                     values_fill = 0.0)
  
  # Jaccard distance and clustering
  jac <- vegdist(mat[-1], method = "jaccard", binary = TRUE)
  jac_mat <- as.matrix(jac)
  colnames(jac_mat) <- mat$year
  rownames(jac_mat) <- mat$year
  
  hc <- hclust(as.dist(jac_mat))
  dend <- as.dendrogram(hc)
  k <- cutree(dend, k = 3, order_clusters_as_data = FALSE)
  
  # Determine which cluster contains early years (1996-2000) - call this "high upwelling"
  early_years <- c("1996", "1997", "1998", "1999", "2000")
  early_cluster <- names(which.max(table(k[early_years])))
  
  # Standardize labeling: Cluster 1 = high upwelling (early years), Cluster 2 = low upwelling
  if (early_cluster == "2") {
    k <- ifelse(k == 1, 2, 1)
  }
  
  list(
    label = label,
    years = mat$year,
    clusters = k,
    dendrogram = dend,
    hclust = hc,
    n_years = nrow(mat),
    n_genera = ncol(mat) - 1
  )
}

# --- 5. RUN ALL THREE ANALYSES ---

# Original: All months, all years
result_original <- run_clustering(
  phyto_counts,
  "Original (all data)"
)

# Option A: High-coverage months only
result_optionA <- run_clustering(
  phyto_counts %>%
    mutate(month = as.numeric(format(date, "%m"))) %>%
    filter(month %in% OPTION_A_MONTHS),
  "Option A (Feb, Mar, Jun, Dec)"
)

# Option B: Well-sampled years only
result_optionB <- run_clustering(
  phyto_counts %>%
    mutate(year = as.numeric(format(date, "%Y"))) %>%
    filter(year %in% OPTION_B_YEARS),
  "Option B (years with ≥10 months)"
)

# --- 6. COMPARE RESULTS ---
cat("\n=== CLUSTERING COMPARISON ===\n\n")

compare_results <- function(res) {
  cat(sprintf("--- %s ---\n", res$label))
  cat(sprintf("Years: %d, Genera: %d\n", res$n_years, res$n_genera))
  
  # Identify cluster membership
  c1_years <- names(res$clusters[res$clusters == 1])
  c2_years <- names(res$clusters[res$clusters == 2])
  
  cat("Cluster 1 (high upwelling):", paste(sort(c1_years), collapse = ", "), "\n")
  cat("Cluster 2 (low upwelling):", paste(sort(c2_years), collapse = ", "), "\n")
  
  # Check key years
  late_years <- c("2014", "2015", "2016")
  late_in_c1 <- sum(late_years %in% c1_years)
  cat(sprintf("Late years (2014-16) in Cluster 1: %d/3\n\n", late_in_c1))
  
  return(data.frame(
    analysis = res$label,
    year = res$years,
    cluster = res$clusters[res$years]
  ))
}

df_original <- compare_results(result_original)
df_optionA <- compare_results(result_optionA)
df_optionB <- compare_results(result_optionB)

# --- 7. SIDE-BY-SIDE COMPARISON TABLE ---
cat("=== YEAR-BY-YEAR COMPARISON ===\n")

comparison <- df_original %>%
  rename(Original = cluster) %>%
  left_join(
    df_optionA %>% rename(OptionA = cluster) %>% select(year, OptionA),
    by = "year"
  ) %>%
  left_join(
    df_optionB %>% rename(OptionB = cluster) %>% select(year, OptionB),
    by = "year"
  ) %>%
  select(year, Original, OptionA, OptionB) %>%
  mutate(
    Consistent = case_when(
      is.na(OptionA) | is.na(OptionB) ~ "—",
      Original == OptionA & Original == OptionB ~ "✓",
      TRUE ~ "✗"
    )
  )

print(as.data.frame(comparison))

# Count consistent assignments
n_consistent <- sum(comparison$Consistent == "✓", na.rm = TRUE)
n_total <- sum(comparison$Consistent != "—", na.rm = TRUE)
cat(sprintf("\nConsistent across all analyses: %d/%d years (%.1f%%)\n", 
            n_consistent, n_total, 100 * n_consistent / n_total))

# --- 8. KEY HYPOTHESIS CHECK ---
cat("\n=== HYPOTHESIS CHECK: Does 2014+ return to Cluster 1? ===\n")

check_return <- function(res) {
  c1 <- names(res$clusters[res$clusters == 1])
  early_in_c1 <- all(c("1996", "1997", "1998", "1999") %in% c1)
  late_in_c1 <- any(c("2014", "2015", "2016") %in% c1)
  middle_in_c2 <- all(c("2005", "2006", "2007", "2008", "2009") %in% 
                        names(res$clusters[res$clusters == 2]))
  
  supports_hypothesis <- early_in_c1 & late_in_c1 & middle_in_c2
  
  cat(sprintf("%s: %s\n", res$label, 
              ifelse(supports_hypothesis, "SUPPORTS hypothesis", "Does NOT support hypothesis")))
  cat(sprintf("  - Early years (96-99) in C1: %s\n", ifelse(early_in_c1, "Yes", "No")))
  cat(sprintf("  - Middle years (05-09) in C2: %s\n", ifelse(middle_in_c2, "Yes", "No")))
  cat(sprintf("  - Any late years (14-16) in C1: %s\n\n", ifelse(late_in_c1, "Yes", "No")))
  
  return(supports_hypothesis)
}

h1 <- check_return(result_original)
h2 <- check_return(result_optionA)
h3 <- check_return(result_optionB)

if (all(c(h1, h2, h3))) {
  cat("==> CONCLUSION: Your hypothesis is ROBUST across sensitivity analyses.\n")
} else if (h1 & (h2 | h3)) {
  cat("==> CONCLUSION: Your hypothesis is PARTIALLY supported. Consider discussing caveats.\n")
} else {
  cat("==> CONCLUSION: Results are SENSITIVE to data filtering. Interpretation requires caution.\n")
}

# --- 9. VISUALIZE DENDROGRAMS SIDE-BY-SIDE ---
# Create color mappings
make_dend_plot <- function(res, title) {
  nm <- setNames(c("Blue", "Red"), c("1", "2"))
  cols <- nm[as.character(res$clusters[res$years])]
  
  dend_styled <- res$dendrogram %>%
    set("labels_cex", 0.8) %>%
    hang.dendrogram(hang_height = 0.1) %>%
    set("leaves_pch", 19) %>%
    set("leaves_col", cols, order_value = TRUE) %>%
    set("branches_lwd", 0.5)
  
  p <- ggplot(dend_styled, horiz = FALSE, offset_labels = -0.03) +
    theme_cowplot(font_size = 12) +
    ylab("Height") +
    ggtitle(title) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 0.75)) +
    scale_x_reverse()
  
  return(p)
}

p1 <- make_dend_plot(result_original, "Original (all data)")
p2 <- make_dend_plot(result_optionA, "Option A (Feb, Mar, Jun, Dec)")
p3 <- make_dend_plot(result_optionB, "Option B (≥10 months/year)")

# Combine plots
combined_plot <- plot_grid(p1, p2, p3, nrow = 1, labels = c("A", "B", "C"))
print(combined_plot)

# Save if desired
# ggsave("plots/Sensitivity_Analysis_Dendrograms.png", combined_plot, width = 15, height = 5, dpi = 300)

cat("\n=== Analysis complete ===\n")

