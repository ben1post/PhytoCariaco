# ============================================
# PRIMARY PRODUCTIVITY ANALYSIS FOR REVIEWER RESPONSE
# ============================================

library(tidyverse)

# Load data (assuming your existing setup)
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")
CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")

# --- 1. DATA COVERAGE ANALYSIS ---
cat("\n=== 1. DATA COVERAGE ===\n")

# Total observations and missing values
n_total <- nrow(CARIACO)
n_pp_available <- sum(!is.na(CARIACO$PrimaryProductivity))
n_pp_missing <- sum(is.na(CARIACO$PrimaryProductivity))
n_chl_available <- sum(!is.na(CARIACO$Chlorophyll))

cat(sprintf("Total monthly observations: %d\n", n_total))
cat(sprintf("PP available: %d (%.1f%%)\n", n_pp_available, 100*n_pp_available/n_total))
cat(sprintf("PP missing: %d (%.1f%%)\n", n_pp_missing, 100*n_pp_missing/n_total))
cat(sprintf("Chlorophyll available: %d (%.1f%%)\n", n_chl_available, 100*n_chl_available/n_total))

# Missing data by year
cat("\n--- PP missing data by year ---\n")
pp_by_year <- CARIACO %>%
  mutate(year = format(date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    n_months = n(),
    pp_available = sum(!is.na(PrimaryProductivity)),
    pp_missing = sum(is.na(PrimaryProductivity)),
    .groups = "drop"
  ) %>%
  filter(pp_missing > 0)
print(pp_by_year, n = 30)

# Count for gradient forest period (matching your n=193 samples)
cat("\n--- Coverage in gradient forest analysis period ---\n")
# Assuming your GF analysis used data where all predictors were available
gf_data <- CARIACO %>%
  filter(!is.na(Chlorophyll) & !is.na(AMO) & !is.na(sst_10m) & !is.na(NO3_merged))
n_gf_total <- nrow(gf_data)
n_gf_with_pp <- sum(!is.na(gf_data$PrimaryProductivity))
cat(sprintf("GF analysis samples (approx): %d\n", n_gf_total))
cat(sprintf("Of these, PP available: %d (%.1f%%)\n", n_gf_with_pp, 100*n_gf_with_pp/n_gf_total))


# --- 2. PP-CHLOROPHYLL CORRELATION ---
cat("\n=== 2. PP-CHLOROPHYLL CORRELATION ===\n")

# Spearman correlation
cor_spearman <- cor.test(CARIACO$PrimaryProductivity, CARIACO$Chlorophyll, 
                         method = "spearman", use = "complete.obs")
cat(sprintf("Spearman rho: %.3f\n", cor_spearman$estimate))
cat(sprintf("Spearman p-value: %.2e\n", cor_spearman$p.value))
cat(sprintf("N pairs: %d\n", sum(!is.na(CARIACO$PrimaryProductivity) & !is.na(CARIACO$Chlorophyll))))

# Pearson on log-transformed values
log_cor <- cor.test(log10(CARIACO$PrimaryProductivity), log10(CARIACO$Chlorophyll), 
                    method = "pearson", use = "complete.obs")
cat(sprintf("\nLog-log Pearson r: %.3f\n", log_cor$estimate))
cat(sprintf("Log-log Pearson p-value: %.2e\n", log_cor$p.value))

# Linear regression log(PP) ~ log(Chl)
lm_fit <- lm(log10(PrimaryProductivity) ~ log10(Chlorophyll), data = CARIACO)
cat(sprintf("\nLinear regression log10(PP) ~ log10(Chl):\n"))
cat(sprintf("  R-squared: %.3f\n", summary(lm_fit)$r.squared))
cat(sprintf("  Slope: %.3f\n", coef(lm_fit)[2]))
cat(sprintf("  Intercept: %.3f\n", coef(lm_fit)[1]))


# --- 3. PP CORRELATIONS WITH KEY PREDICTORS ---
cat("\n=== 3. PP CORRELATIONS WITH KEY PREDICTORS ===\n")

predictors <- c("AMO", "sst_10m", "NO3_merged", "MEIv2", "Isotherm_21")
predictor_names <- c("AMO", "SST", "NO3", "MEI v.2", "21°C Isotherm")

for (i in seq_along(predictors)) {
  cor_test <- cor.test(CARIACO$PrimaryProductivity, CARIACO[[predictors[i]]], 
                       method = "spearman", use = "complete.obs")
  cat(sprintf("%s: rho = %.3f, p = %.2e\n", 
              predictor_names[i], cor_test$estimate, cor_test$p.value))
}


# --- 4. WILCOXON TEST: PP BETWEEN CLUSTERS ---
cat("\n=== 4. WILCOXON TEST: PP BETWEEN CLUSTERS ===\n")

# Use your existing cluster assignment function
extractGroups <- function(year1, year2, year3, year4){
  Group1 <- CARIACO %>%  
    filter(date >= as.Date(year1, format="%Y-%m-%d") & date <= as.Date(year2, format="%Y-%m-%d"))
  Group2 <- CARIACO %>%  
    filter(date >= as.Date(year3, format="%Y-%m-%d") & date <= as.Date(year4, format="%Y-%m-%d"))
  Group = rbind(Group1, Group2)
  return(Group)
}

GroupRed = extractGroups(year1="1996-01-01", year2="2003-12-31", year3="2014-01-01", year4="2016-12-31")
GroupBlue = extractGroups(year1="2004-01-01", year2="2013-12-31", year3="2017-06-01", year4="2017-12-31")
GroupRed$group = "Cluster 1"
GroupBlue$group = "Cluster 2"
ENV_DATA_groups = rbind(GroupRed, GroupBlue)

# Wilcoxon test for PP
pp_wilcox <- wilcox.test(PrimaryProductivity ~ group, data = ENV_DATA_groups)
cat(sprintf("W statistic: %.1f\n", pp_wilcox$statistic))
cat(sprintf("p-value: %.4f\n", pp_wilcox$p.value))

# Summary statistics by cluster
cat("\n--- PP summary by cluster ---\n")
pp_summary <- ENV_DATA_groups %>%
  group_by(group) %>%
  summarise(
    n = sum(!is.na(PrimaryProductivity)),
    n_missing = sum(is.na(PrimaryProductivity)),
    mean_PP = mean(PrimaryProductivity, na.rm = TRUE),
    median_PP = median(PrimaryProductivity, na.rm = TRUE),
    sd_PP = sd(PrimaryProductivity, na.rm = TRUE),
    .groups = "drop"
  )
print(pp_summary)

# For comparison: Chlorophyll Wilcoxon (should match your Table 2)
chl_wilcox <- wilcox.test(Chlorophyll ~ group, data = ENV_DATA_groups)
cat(sprintf("\nFor comparison - Chlorophyll a:\n"))
cat(sprintf("W statistic: %.1f, p-value: %.2e\n", chl_wilcox$statistic, chl_wilcox$p.value))


# --- 5. COMPARISON TABLE ---
cat("\n=== 5. SUMMARY FOR RESPONSE ===\n")
cat("Copy these values into your reviewer response:\n\n")
cat(sprintf("- PP-Chl Spearman correlation: ρ = %.2f, p < 0.001\n", cor_spearman$estimate))
cat(sprintf("- Log-log R² = %.2f (%.0f%% variance explained)\n", 
            summary(lm_fit)$r.squared, 100*summary(lm_fit)$r.squared))
cat(sprintf("- PP data coverage: %d of %d months (%.0f%% missing)\n", 
            n_pp_available, n_total, 100*n_pp_missing/n_total))
cat(sprintf("- Wilcoxon test PP between clusters: W = %.0f, p = %.3f\n", 
            pp_wilcox$statistic, pp_wilcox$p.value))
