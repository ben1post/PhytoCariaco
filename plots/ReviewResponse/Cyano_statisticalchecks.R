# ============================================
# CYANOBACTERIA SEASONAL/CLUSTER ANALYSIS
# ============================================

library(tidyverse)

# Load data and prep (using your existing setup)
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  factor(season_lookup[months], levels = c("Upwelling", "Secondary Upwelling", "Rainy"))
}

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    season = get_season_3cat(date),
    cluster = case_when(
      (date >= as.Date("1996-01-01") & date <= as.Date("2003-12-31")) |
        (date >= as.Date("2014-01-01") & date <= as.Date("2016-12-31")) ~ "Cluster 1",
      (date >= as.Date("2004-01-01") & date <= as.Date("2013-12-31")) |
        (date >= as.Date("2017-06-01") & date <= as.Date("2017-12-31")) ~ "Cluster 2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cluster))


# --- 1. DESCRIPTIVE STATISTICS ---
cat("\n=== 1. DESCRIPTIVE STATISTICS: Cyanobacteria by Season and Cluster ===\n")

cyano_summary <- CARIACO %>%
  group_by(cluster, season) %>%
  summarise(
    n = sum(!is.na(Abundance_Cyano)),
    mean = mean(Abundance_Cyano, na.rm = TRUE),
    median = median(Abundance_Cyano, na.rm = TRUE),
    sd = sd(Abundance_Cyano, na.rm = TRUE),
    .groups = "drop"
  )
print(cyano_summary)


# --- 2. SEASONAL DIFFERENCES WITHIN EACH CLUSTER ---
cat("\n=== 2. KRUSKAL-WALLIS: Seasonal differences within each cluster ===\n")

# Cluster 1
kw_c1 <- kruskal.test(Abundance_Cyano ~ season, data = CARIACO %>% filter(cluster == "Cluster 1"))
cat(sprintf("Cluster 1 - Kruskal-Wallis: chi-squared = %.2f, df = %d, p = %.4f\n", 
            kw_c1$statistic, kw_c1$parameter, kw_c1$p.value))

# Cluster 2
kw_c2 <- kruskal.test(Abundance_Cyano ~ season, data = CARIACO %>% filter(cluster == "Cluster 2"))
cat(sprintf("Cluster 2 - Kruskal-Wallis: chi-squared = %.2f, df = %d, p = %.4f\n", 
            kw_c2$statistic, kw_c2$parameter, kw_c2$p.value))


# --- 3. CLUSTER DIFFERENCES WITHIN EACH SEASON ---
cat("\n=== 3. WILCOXON: Cluster differences within each season ===\n")

seasons <- c("Upwelling", "Secondary Upwelling", "Rainy")
for (s in seasons) {
  test_data <- CARIACO %>% filter(season == s)
  wt <- wilcox.test(Abundance_Cyano ~ cluster, data = test_data)
  
  # Get medians for context
  med_c1 <- median(test_data$Abundance_Cyano[test_data$cluster == "Cluster 1"], na.rm = TRUE)
  med_c2 <- median(test_data$Abundance_Cyano[test_data$cluster == "Cluster 2"], na.rm = TRUE)
  
  cat(sprintf("%s: W = %.0f, p = %.4f (Median C1: %.2f, C2: %.2f)\n", 
              s, wt$statistic, wt$p.value, med_c1, med_c2))
}


# --- 4. INTERACTION TEST (Scheirer-Ray-Hare) ---
cat("\n=== 4. SCHEIRER-RAY-HARE TEST: Cluster × Season interaction ===\n")

# This is a non-parametric two-way ANOVA equivalent
# Using the rcompanion package if available, otherwise manual calculation
if (requireNamespace("rcompanion", quietly = TRUE)) {
  library(rcompanion)
  srh <- scheirerRayHare(Abundance_Cyano ~ cluster + season, data = CARIACO)
  print(srh)
} else {
  cat("Note: Install 'rcompanion' package for Scheirer-Ray-Hare test\n")
  cat("Running manual rank-based analysis instead...\n\n")
  
  # Manual approach using ranks
  CARIACO$rank_cyano <- rank(CARIACO$Abundance_Cyano, na.last = "keep")
  aov_ranks <- aov(rank_cyano ~ cluster * season, data = CARIACO)
  cat("ANOVA on ranks (approximation):\n")
  print(summary(aov_ranks))
}


# --- 5. GENUS-LEVEL ANALYSIS (if available) ---
cat("\n=== 5. GENUS-LEVEL: Trichodesmium and Synechococcus ===\n")

# Check if these columns exist in your data
cyano_genera <- c("Trichodesmium", "Synechococcus")

# First, let's see what cyanobacteria-related columns might exist
cyano_cols <- names(CARIACO)[grep("Tricho|Synecho|cyano|Cyano", names(CARIACO), ignore.case = TRUE)]
cat("Available cyanobacteria-related columns:\n")
print(cyano_cols)

# If genus-level data exists in your microscopy data, you may need to load it separately
# For now, report what's available


# --- 6. CORRELATION WITH STRATIFICATION INDICATORS ---
cat("\n=== 6. SPEARMAN CORRELATIONS: Cyanobacteria vs stratification indicators ===\n")

strat_vars <- c("sst_10m", "Isotherm_21", "NO3_merged", "AMO")
strat_names <- c("SST", "21°C Isotherm", "NO3", "AMO")

for (i in seq_along(strat_vars)) {
  cor_test <- cor.test(CARIACO$Abundance_Cyano, CARIACO[[strat_vars[i]]], 
                       method = "spearman", use = "complete.obs")
  cat(sprintf("%s: rho = %.3f, p = %.4f\n", 
              strat_names[i], cor_test$estimate, cor_test$p.value))
}


# --- 7. SUMMARY OUTPUT ---
cat("\n=== 7. SUMMARY FOR REVIEWER RESPONSE ===\n")

# Overall cluster effect (matching your Table 2 approach)
overall_wilcox <- wilcox.test(Abundance_Cyano ~ cluster, data = CARIACO)
cat(sprintf("Overall cluster effect: W = %.0f, p = %.4f\n", 
            overall_wilcox$statistic, overall_wilcox$p.value))

cat("\nKey findings to report:\n")
cat("- Seasonal effect within Cluster 1: p = ", sprintf("%.4f", kw_c1$p.value), "\n")
cat("- Seasonal effect within Cluster 2: p = ", sprintf("%.4f", kw_c2$p.value), "\n")
cat("- Cluster effect is [stronger/similar/weaker] during rainy season vs upwelling\n")


########## CHECK TRICHODESMIUM AND SYNECHOCOCCUS INDIVIDUALLY ##############


# ============================================
# GENUS-LEVEL CYANOBACTERIA ANALYSIS
# Trichodesmium and Synechococcus
# ============================================

library(tidyverse)

# --- 1. LOAD DATA ---
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts),]

CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# --- 2. PREP CARIACO WITH CLUSTERS AND SEASONS ---
get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  factor(season_lookup[months], levels = c("Upwelling", "Secondary Upwelling", "Rainy"))
}

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    season = get_season_3cat(date),
    cluster = case_when(
      (date >= as.Date("1996-01-01") & date <= as.Date("2003-12-31")) |
        (date >= as.Date("2014-01-01") & date <= as.Date("2016-12-31")) ~ "Cluster 1",
      (date >= as.Date("2004-01-01") & date <= as.Date("2013-12-31")) |
        (date >= as.Date("2017-06-01") & date <= as.Date("2017-12-31")) ~ "Cluster 2",
      TRUE ~ NA_character_
    )
  )

# --- 3. EXTRACT AND AGGREGATE TRICHODESMIUM & SYNECHOCOCCUS ---
cyano_genera <- phyto_counts %>%
  filter(Genus %in% c("Trichodesmium", "Synechococcus")) %>%
  group_by(date, Genus) %>%
  summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = counts) %>%
  # Convert to year-month Date (15th of each month) to match CARIACO
  
  mutate(date = as.Date(format(date, "%Y-%m-15")))

# --- 4. MERGE WITH CARIACO ---
cyano_data <- CARIACO %>%
  select(date, cluster, season, sst_10m, Isotherm_21, NO3_merged, AMO) %>%
  inner_join(cyano_genera, by = "date") %>%
  filter(!is.na(cluster))

cat("\n=== DATA OVERVIEW ===\n")
cat(sprintf("Total observations with phyto data: %d\n", nrow(cyano_data)))
cat(sprintf("Trichodesmium non-NA: %d\n", sum(!is.na(cyano_data$Trichodesmium))))
cat(sprintf("Synechococcus non-NA: %d\n", sum(!is.na(cyano_data$Synechococcus))))


# --- 5. DESCRIPTIVE STATISTICS ---
cat("\n=== DESCRIPTIVE STATISTICS ===\n")

for (genus in c("Trichodesmium", "Synechococcus")) {
  cat(sprintf("\n--- %s ---\n", genus))
  
  summary_tbl <- cyano_data %>%
    group_by(cluster, season) %>%
    summarise(
      n_obs = sum(!is.na(.data[[genus]])),
      mean = mean(.data[[genus]], na.rm = TRUE),
      median = median(.data[[genus]], na.rm = TRUE),
      max = max(.data[[genus]], na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_tbl)
}


# --- 6. STATISTICAL TESTS: TRICHODESMIUM ---
cat("\n\n========== TRICHODESMIUM ==========\n")

# Filter to non-NA observations for this genus
tricho_data <- cyano_data %>% filter(!is.na(Trichodesmium))
cat(sprintf("N observations: %d\n", nrow(tricho_data)))

# Overall cluster effect
wt_tricho_overall <- wilcox.test(Trichodesmium ~ cluster, data = tricho_data)
cat(sprintf("\nOverall cluster effect: W = %.0f, p = %.4f\n", 
            wt_tricho_overall$statistic, wt_tricho_overall$p.value))

# Seasonal differences within clusters
cat("\nKruskal-Wallis (seasonal effect within clusters):\n")
tricho_c1 <- tricho_data %>% filter(cluster == "Cluster 1")
tricho_c2 <- tricho_data %>% filter(cluster == "Cluster 2")

if (nrow(tricho_c1) >= 3) {
  kw_tricho_c1 <- kruskal.test(Trichodesmium ~ season, data = tricho_c1)
  cat(sprintf("  Cluster 1 (n=%d): chi-sq = %.2f, p = %.4f\n", 
              nrow(tricho_c1), kw_tricho_c1$statistic, kw_tricho_c1$p.value))
} else {
  cat(sprintf("  Cluster 1: insufficient data (n=%d)\n", nrow(tricho_c1)))
}

if (nrow(tricho_c2) >= 3) {
  kw_tricho_c2 <- kruskal.test(Trichodesmium ~ season, data = tricho_c2)
  cat(sprintf("  Cluster 2 (n=%d): chi-sq = %.2f, p = %.4f\n", 
              nrow(tricho_c2), kw_tricho_c2$statistic, kw_tricho_c2$p.value))
} else {
  cat(sprintf("  Cluster 2: insufficient data (n=%d)\n", nrow(tricho_c2)))
}

# Cluster differences within seasons
cat("\nWilcoxon (cluster effect within seasons):\n")
for (s in c("Upwelling", "Secondary Upwelling", "Rainy")) {
  test_data <- tricho_data %>% filter(season == s)
  n_c1 <- sum(test_data$cluster == "Cluster 1")
  n_c2 <- sum(test_data$cluster == "Cluster 2")
  
  if (n_c1 >= 1 & n_c2 >= 1) {
    wt <- wilcox.test(Trichodesmium ~ cluster, data = test_data)
    med_c1 <- median(test_data$Trichodesmium[test_data$cluster == "Cluster 1"], na.rm = TRUE)
    med_c2 <- median(test_data$Trichodesmium[test_data$cluster == "Cluster 2"], na.rm = TRUE)
    cat(sprintf("  %s (n=%d/%d): W = %.0f, p = %.4f (Median C1: %.3f, C2: %.3f)\n", 
                s, n_c1, n_c2, wt$statistic, wt$p.value, med_c1, med_c2))
  } else {
    cat(sprintf("  %s: insufficient data (n=%d/%d)\n", s, n_c1, n_c2))
  }
}


# --- 7. STATISTICAL TESTS: SYNECHOCOCCUS ---
cat("\n\n========== SYNECHOCOCCUS ==========\n")

# Filter to non-NA observations for this genus
synecho_data <- cyano_data %>% filter(!is.na(Synechococcus))
cat(sprintf("N observations: %d\n", nrow(synecho_data)))

# Overall cluster effect
wt_synecho_overall <- wilcox.test(Synechococcus ~ cluster, data = synecho_data)
cat(sprintf("\nOverall cluster effect: W = %.0f, p = %.4f\n", 
            wt_synecho_overall$statistic, wt_synecho_overall$p.value))

# Seasonal differences within clusters
cat("\nKruskal-Wallis (seasonal effect within clusters):\n")
synecho_c1 <- synecho_data %>% filter(cluster == "Cluster 1")
synecho_c2 <- synecho_data %>% filter(cluster == "Cluster 2")

if (nrow(synecho_c1) >= 3) {
  kw_synecho_c1 <- kruskal.test(Synechococcus ~ season, data = synecho_c1)
  cat(sprintf("  Cluster 1 (n=%d): chi-sq = %.2f, p = %.4f\n", 
              nrow(synecho_c1), kw_synecho_c1$statistic, kw_synecho_c1$p.value))
} else {
  cat(sprintf("  Cluster 1: insufficient data (n=%d)\n", nrow(synecho_c1)))
}

if (nrow(synecho_c2) >= 3) {
  kw_synecho_c2 <- kruskal.test(Synechococcus ~ season, data = synecho_c2)
  cat(sprintf("  Cluster 2 (n=%d): chi-sq = %.2f, p = %.4f\n", 
              nrow(synecho_c2), kw_synecho_c2$statistic, kw_synecho_c2$p.value))
} else {
  cat(sprintf("  Cluster 2: insufficient data (n=%d)\n", nrow(synecho_c2)))
}

# Cluster differences within seasons
cat("\nWilcoxon (cluster effect within seasons):\n")
for (s in c("Upwelling", "Secondary Upwelling", "Rainy")) {
  test_data <- synecho_data %>% filter(season == s)
  n_c1 <- sum(test_data$cluster == "Cluster 1")
  n_c2 <- sum(test_data$cluster == "Cluster 2")
  
  if (n_c1 >= 1 & n_c2 >= 1) {
    wt <- wilcox.test(Synechococcus ~ cluster, data = test_data)
    med_c1 <- median(test_data$Synechococcus[test_data$cluster == "Cluster 1"], na.rm = TRUE)
    med_c2 <- median(test_data$Synechococcus[test_data$cluster == "Cluster 2"], na.rm = TRUE)
    cat(sprintf("  %s (n=%d/%d): W = %.0f, p = %.4f (Median C1: %.3f, C2: %.3f)\n", 
                s, n_c1, n_c2, wt$statistic, wt$p.value, med_c1, med_c2))
  } else {
    cat(sprintf("  %s: insufficient data (n=%d/%d)\n", s, n_c1, n_c2))
  }
}


# --- 8. CORRELATIONS WITH STRATIFICATION ---
cat("\n\n=== CORRELATIONS WITH STRATIFICATION INDICATORS ===\n")

strat_vars <- c("sst_10m", "Isotherm_21", "NO3_merged", "AMO")
strat_names <- c("SST", "21°C Isotherm", "NO3", "AMO")

cat("\nTrichodesmium:\n")
for (i in seq_along(strat_vars)) {
  ct <- cor.test(tricho_data$Trichodesmium, tricho_data[[strat_vars[i]]], 
                 method = "spearman", use = "complete.obs")
  n_pairs <- sum(complete.cases(tricho_data$Trichodesmium, tricho_data[[strat_vars[i]]]))
  cat(sprintf("  %s (n=%d): rho = %.3f, p = %.4f\n", strat_names[i], n_pairs, ct$estimate, ct$p.value))
}

cat("\nSynechococcus:\n")
for (i in seq_along(strat_vars)) {
  ct <- cor.test(synecho_data$Synechococcus, synecho_data[[strat_vars[i]]], 
                 method = "spearman", use = "complete.obs")
  n_pairs <- sum(complete.cases(synecho_data$Synechococcus, synecho_data[[strat_vars[i]]]))
  cat(sprintf("  %s (n=%d): rho = %.3f, p = %.4f\n", strat_names[i], n_pairs, ct$estimate, ct$p.value))
}


# --- 9. SUMMARY ---
cat("\n\n=== SUMMARY FOR REVIEWER RESPONSE ===\n")

cat("\nTrichodesmium:\n")
cat(sprintf("  - N observations: %d\n", nrow(tricho_data)))
cat(sprintf("  - Overall cluster effect: p = %.4f\n", wt_tricho_overall$p.value))
med_tricho_c1 <- median(tricho_data$Trichodesmium[tricho_data$cluster == "Cluster 1"], na.rm = TRUE)
med_tricho_c2 <- median(tricho_data$Trichodesmium[tricho_data$cluster == "Cluster 2"], na.rm = TRUE)
cat(sprintf("  - Median Cluster 1: %.4f, Cluster 2: %.4f\n", med_tricho_c1, med_tricho_c2))

cat("\nSynechococcus:\n")
cat(sprintf("  - N observations: %d\n", nrow(synecho_data)))
cat(sprintf("  - Overall cluster effect: p = %.4f\n", wt_synecho_overall$p.value))
med_synecho_c1 <- median(synecho_data$Synechococcus[synecho_data$cluster == "Cluster 1"], na.rm = TRUE)
med_synecho_c2 <- median(synecho_data$Synechococcus[synecho_data$cluster == "Cluster 2"], na.rm = TRUE)
cat(sprintf("  - Median Cluster 1: %.4f, Cluster 2: %.4f\n", med_synecho_c1, med_synecho_c2))

