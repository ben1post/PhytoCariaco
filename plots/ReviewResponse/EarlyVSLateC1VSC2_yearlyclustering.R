# ============================================
# PART 2: YEAR-BY-YEAR SIMILARITY ANALYSIS
# ============================================

library(tidyverse)
library(vegan)

# --- SWITCH: Core months or all months ---
USE_CORE_MONTHS <- TRUE  # Set to TRUE to restrict to core months

CORE_MONTHS <- c(1, 2, 4, 5, 9, 11, 12)

# --- 1. LOAD DATA ---
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# --- 2. PREP DATA ---
phyto_data <- phyto_counts %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  mutate(
    year = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m"))
  ) %>%
  filter(year >= 1996 & year <= 2016)  # Match original clustering years

# --- 3. APPLY CORE MONTHS FILTER (if enabled) ---
if (USE_CORE_MONTHS) {
  phyto_data <- phyto_data %>% filter(month %in% CORE_MONTHS)
  cat("=== USING CORE MONTHS ONLY ===\n\n")
} else {
  cat("=== USING ALL AVAILABLE MONTHS ===\n\n")
}

# --- 4. AGGREGATE BY YEAR ---
year_genus <- phyto_data %>%
  group_by(year, Genus) %>%
  summarise(Total = sum(counts, na.rm = TRUE), .groups = "drop")

year_matrix <- year_genus %>%
  pivot_wider(names_from = Genus, values_from = Total, values_fill = 0) %>%
  column_to_rownames("year")

cat("=== YEAR MATRIX DIMENSIONS ===\n")
cat(sprintf("Years: %d, Genera: %d\n\n", nrow(year_matrix), ncol(year_matrix)))

# --- 5. DEFINE PERIODS ---
early_c1_years <- as.character(1996:2004)
c2_years <- as.character(2005:2013)
late_c1_years <- as.character(2014:2016)

cat("=== PERIODS ===\n")
cat("Early C1:", paste(early_c1_years, collapse = ", "), "\n")
cat("C2:", paste(c2_years, collapse = ", "), "\n")
cat("Late C1:", paste(late_c1_years, collapse = ", "), "\n\n")

# --- 6. CALCULATE DISTANCE MATRICES ---
jaccard_dist <- vegdist(year_matrix, method = "jaccard", binary = TRUE)
bray_dist <- vegdist(year_matrix, method = "bray", binary = FALSE)

# Convert to similarity matrices
jaccard_sim <- 1 - as.matrix(jaccard_dist)
bray_sim <- 1 - as.matrix(bray_dist)

# --- 7. HELPER FUNCTION: Extract pairwise similarities ---
get_pairwise_similarities <- function(sim_matrix, years1, years2 = NULL) {
  # If years2 is NULL, get within-group similarities
  if (is.null(years2)) {
    years2 <- years1
    within_group <- TRUE
  } else {
    within_group <- FALSE
  }
  
  # Filter to available years
  years1 <- years1[years1 %in% rownames(sim_matrix)]
  years2 <- years2[years2 %in% rownames(sim_matrix)]
  
  # Extract submatrix
  sub_matrix <- sim_matrix[years1, years2, drop = FALSE]
  
  if (within_group) {
    # Get upper triangle only (exclude diagonal and duplicates)
    values <- sub_matrix[upper.tri(sub_matrix)]
  } else {
    # Get all pairwise values
    values <- as.vector(sub_matrix)
  }
  
  return(values)
}

# --- 8. CALCULATE MEAN SIMILARITIES ---
cat("=== JACCARD SIMILARITY (Presence/Absence) ===\n")

j_within_early <- get_pairwise_similarities(jaccard_sim, early_c1_years)
j_within_c2 <- get_pairwise_similarities(jaccard_sim, c2_years)
j_within_late <- get_pairwise_similarities(jaccard_sim, late_c1_years)
j_late_vs_early <- get_pairwise_similarities(jaccard_sim, late_c1_years, early_c1_years)
j_late_vs_c2 <- get_pairwise_similarities(jaccard_sim, late_c1_years, c2_years)

cat(sprintf("Within Early C1:      %.3f ± %.3f (n=%d pairs)\n", 
            mean(j_within_early), sd(j_within_early), length(j_within_early)))
cat(sprintf("Within C2:            %.3f ± %.3f (n=%d pairs)\n", 
            mean(j_within_c2), sd(j_within_c2), length(j_within_c2)))
cat(sprintf("Within Late C1:       %.3f ± %.3f (n=%d pairs)\n", 
            mean(j_within_late), sd(j_within_late), length(j_within_late)))
cat(sprintf("Late C1 vs Early C1:  %.3f ± %.3f (n=%d pairs)\n", 
            mean(j_late_vs_early), sd(j_late_vs_early), length(j_late_vs_early)))
cat(sprintf("Late C1 vs C2:        %.3f ± %.3f (n=%d pairs)\n", 
            mean(j_late_vs_c2), sd(j_late_vs_c2), length(j_late_vs_c2)))

cat("\n=== BRAY-CURTIS SIMILARITY (Abundance) ===\n")

b_within_early <- get_pairwise_similarities(bray_sim, early_c1_years)
b_within_c2 <- get_pairwise_similarities(bray_sim, c2_years)
b_within_late <- get_pairwise_similarities(bray_sim, late_c1_years)
b_late_vs_early <- get_pairwise_similarities(bray_sim, late_c1_years, early_c1_years)
b_late_vs_c2 <- get_pairwise_similarities(bray_sim, late_c1_years, c2_years)

cat(sprintf("Within Early C1:      %.3f ± %.3f (n=%d pairs)\n", 
            mean(b_within_early), sd(b_within_early), length(b_within_early)))
cat(sprintf("Within C2:            %.3f ± %.3f (n=%d pairs)\n", 
            mean(b_within_c2), sd(b_within_c2), length(b_within_c2)))
cat(sprintf("Within Late C1:       %.3f ± %.3f (n=%d pairs)\n", 
            mean(b_within_late), sd(b_within_late), length(b_within_late)))
cat(sprintf("Late C1 vs Early C1:  %.3f ± %.3f (n=%d pairs)\n", 
            mean(b_late_vs_early), sd(b_late_vs_early), length(b_late_vs_early)))
cat(sprintf("Late C1 vs C2:        %.3f ± %.3f (n=%d pairs)\n", 
            mean(b_late_vs_c2), sd(b_late_vs_c2), length(b_late_vs_c2)))

# --- 9. STATISTICAL TESTS ---
cat("\n=== STATISTICAL TESTS ===\n")
cat("Testing: Is Late C1 more similar to Early C1 or to C2?\n\n")

# Wilcoxon test comparing Late-vs-Early similarities to Late-vs-C2 similarities
cat("Jaccard:\n")
wt_jaccard <- wilcox.test(j_late_vs_early, j_late_vs_c2)
cat(sprintf("  Late vs Early C1 mean: %.3f\n", mean(j_late_vs_early)))
cat(sprintf("  Late vs C2 mean:       %.3f\n", mean(j_late_vs_c2)))
cat(sprintf("  Wilcoxon p-value:      %.4f\n", wt_jaccard$p.value))
if (mean(j_late_vs_early) > mean(j_late_vs_c2)) {
  cat("  → Late C1 is MORE similar to Early C1 (supports return)\n")
} else {
  cat("  → Late C1 is MORE similar to C2 (does not support return)\n")
}

cat("\nBray-Curtis:\n")
wt_bray <- wilcox.test(b_late_vs_early, b_late_vs_c2)
cat(sprintf("  Late vs Early C1 mean: %.3f\n", mean(b_late_vs_early)))
cat(sprintf("  Late vs C2 mean:       %.3f\n", mean(b_late_vs_c2)))
cat(sprintf("  Wilcoxon p-value:      %.4f\n", wt_bray$p.value))
if (mean(b_late_vs_early) > mean(b_late_vs_c2)) {
  cat("  → Late C1 is MORE similar to Early C1 (supports return)\n")
} else {
  cat("  → Late C1 is MORE similar to C2 (does not support return)\n")
}

# --- 10. SUMMARY ---
cat("\n=== SUMMARY ===\n")
cat("Composition (Jaccard):\n")
if (mean(j_late_vs_early) > mean(j_late_vs_c2)) {
  cat("  ✓ Late C1 years are on average MORE similar to Early C1 years\n")
} else {
  cat("  ✗ Late C1 years are on average MORE similar to C2 years\n")
}

cat("Dominance (Bray-Curtis):\n
")
if (mean(b_late_vs_early) > mean(b_late_vs_c2)) {
  cat("  ✓ Late C1 years are on average MORE similar to Early C1 years\n")
} else {
  cat("  ✗ Late C1 years are on average MORE similar to C2 years\n")
}