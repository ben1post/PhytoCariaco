# ============================================
# CORE MONTHS SENSITIVITY ANALYSIS
# ============================================

library(tidyverse)

# --- 1. LOAD AND PREP DATA ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    year = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m")),
    period = case_when(
      year >= 1996 & year <= 2004 ~ "Early_C1",
      year >= 2005 & year <= 2013 ~ "C2",
      year >= 2014 & year <= 2017 ~ "Late_C1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period))


# --- 2. IDENTIFY CORE MONTHS ---
# Find months with phytoplankton data in ALL THREE periods

# Using Chlorophyll as proxy for data availability (similar pattern to phyto counts)
months_by_period <- CARIACO %>%
  filter(!is.na(Chlorophyll)) %>%
  group_by(period, month) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = n_obs, values_fill = 0)

cat("=== MONTHS AVAILABLE BY PERIOD ===\n")
print(months_by_period)

# Core months: present in all three periods
core_months <- months_by_period %>%
  filter(Early_C1 > 0 & C2 > 0 & Late_C1 > 0) %>%
  pull(month)

cat("\nCore months available in all periods:", paste(month.abb[core_months], collapse = ", "), "\n")
cat("Number of core months:", length(core_months), "\n")


# --- 3. FILTER TO CORE MONTHS ---
CARIACO_core <- CARIACO %>%
  filter(month %in% core_months)

cat("\n=== SAMPLE SIZES (Core Months Only) ===\n")
CARIACO_core %>%
  group_by(period) %>%
  summarise(
    n_total = n(),
    n_chl = sum(!is.na(Chlorophyll)),
    n_sst = sum(!is.na(sst_10m)),
    n_no3 = sum(!is.na(NO3_merged)),
    years = paste(range(year), collapse = "-"),
    .groups = "drop"
  ) %>%
  print()


# --- 4. STATISTICAL TESTS ---
cat("\n=== WILCOXON TESTS: CORE MONTHS ONLY ===\n")

# Variables to test
test_vars <- c("sst_10m", "NO3_merged", "Chlorophyll", "Abundance_Cyano", 
               "Abundance_Diatom", "Abundance_Hapto", "GenusRichness")
var_names <- c("SST", "NO3", "Chlorophyll", "Cyanobacteria", 
               "Diatoms", "Haptophytes", "Genus Richness")

# Comparisons
comparisons <- list(
  c("Early_C1", "C2"),
  c("C2", "Late_C1"),
  c("Early_C1", "Late_C1")
)
comp_names <- c("Early_C1 vs C2", "C2 vs Late_C1", "Early_C1 vs Late_C1")

# Run tests
results <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  vname <- var_names[i]
  
  for (j in seq_along(comparisons)) {
    comp <- comparisons[[j]]
    cname <- comp_names[j]
    
    test_data <- CARIACO_core %>%
      filter(period %in% comp) %>%
      filter(!is.na(.data[[var]]))
    
    if (nrow(test_data) > 0 & length(unique(test_data$period)) == 2) {
      # Get sample sizes and medians
      n1 <- sum(test_data$period == comp[1])
      n2 <- sum(test_data$period == comp[2])
      med1 <- median(test_data[[var]][test_data$period == comp[1]], na.rm = TRUE)
      med2 <- median(test_data[[var]][test_data$period == comp[2]], na.rm = TRUE)
      
      # Wilcoxon test
      wt <- wilcox.test(as.formula(paste(var, "~ period")), data = test_data)
      
      results <- rbind(results, data.frame(
        Variable = vname,
        Comparison = cname,
        n1 = n1,
        n2 = n2,
        Median1 = round(med1, 4),
        Median2 = round(med2, 4),
        W = wt$statistic,
        p_value = wt$p.value,
        Sig = case_when(
          wt$p.value < 0.001 ~ "***",
          wt$p.value < 0.01 ~ "**",
          wt$p.value < 0.05 ~ "*",
          TRUE ~ "n.s."
        )
      ))
    }
  }
}

# Print results
cat("\n")
print(results, row.names = FALSE)


# --- 5. SUMMARY TABLE BY PERIOD ---
cat("\n\n=== MEDIANS BY PERIOD (Core Months Only) ===\n")

summary_table <- CARIACO_core %>%
  group_by(period) %>%
  summarise(
    n = n(),
    SST = median(sst_10m, na.rm = TRUE),
    NO3 = median(NO3_merged, na.rm = TRUE),
    Chl = median(Chlorophyll, na.rm = TRUE),
    Diatoms = median(Abundance_Diatom, na.rm = TRUE),
    Cyano = median(Abundance_Cyano, na.rm = TRUE),
    Richness = median(GenusRichness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Reorder periods
  
  mutate(period = factor(period, levels = c("Early_C1", "C2", "Late_C1"))) %>%
  arrange(period)

print(summary_table)


# --- 6. INTERPRETATION HELPER ---
cat("\n\n=== INTERPRETATION SUMMARY ===\n")
cat("If 'return' hypothesis is supported:\n")
cat("  - Early_C1 vs C2: significant differences (original regime shift)\n")
cat("  - C2 vs Late_C1: significant differences (late ≠ low-upwelling)\n")
cat("  - Early_C1 vs Late_C1: non-significant or smaller differences\n")
