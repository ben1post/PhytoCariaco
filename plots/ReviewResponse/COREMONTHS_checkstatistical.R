# ============================================
# DIAGNOSTIC: Monthly Sampling Coverage
# ============================================

library(tidyverse)

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

# --- 3. VISUALIZE: Heatmap of sampling coverage ---
coverage_matrix <- monthly_coverage %>%
  mutate(sampled = as.numeric(sampled)) %>%
  pivot_wider(names_from = month, values_from = sampled) %>%
  column_to_rownames("year")

cat("=== MONTHLY SAMPLING MATRIX (1 = sampled, 0 = missing) ===\n")
print(as.matrix(coverage_matrix))

# Heatmap visualization
p_heatmap <- ggplot(monthly_coverage, aes(x = factor(month), y = factor(year), fill = sampled)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "steelblue"), 
                    labels = c("Missing", "Sampled")) +
  labs(x = "Month", y = "Year", fill = "Status",
       title = "Monthly Sampling Coverage (1996-2016)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())
print(p_heatmap)

# --- 4. IDENTIFY CONSISTENTLY SAMPLED MONTHS ---
month_consistency <- monthly_coverage %>%
  group_by(month) %>%
  summarise(
    years_sampled = sum(sampled),
    total_years = n(),
    pct_coverage = round(100 * years_sampled / total_years, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_coverage))

cat("\n=== MONTH-LEVEL CONSISTENCY (across all years) ===\n")
print(month_consistency, n = 12)

# --- 5. FOCUS ON FINAL YEARS (where you suspect gaps) ---
final_years <- c(2014, 2015, 2016)

final_year_coverage <- monthly_coverage %>%
  filter(year %in% final_years) %>%
  group_by(month) %>%
  summarise(
    years_sampled = sum(sampled),
    years_checked = n(),
    .groups = "drop"
  ) %>%
  mutate(
    all_final_years_sampled = years_sampled == length(final_years),
    month_name = month.abb[month]
  )

cat("\n=== COVERAGE IN FINAL YEARS (2014-2016) ===\n")
print(final_year_coverage)

# Months present in ALL final years:
months_in_all_final <- final_year_coverage %>%
  filter(all_final_years_sampled) %>%
  pull(month)

cat("\n=== MONTHS SAMPLED IN ALL FINAL YEARS ===\n")
cat(paste(month.abb[months_in_all_final], collapse = ", "), "\n")

# --- 6. CHECK: Does using these months create gaps in OTHER years? ---
cat("\n=== CHECKING PROPOSED CORE MONTHS ===\n")

# Your current selection:
CORE_MONTHS_CURRENT <- c(1, 2, 4, 5, 9, 11, 12)

# Data-driven selection (months in all final years):
CORE_MONTHS_DATADRIVEN <- months_in_all_final

check_core_months <- function(core_months, label) {
  cat(sprintf("\n--- %s: months %s ---\n", label, paste(core_months, collapse = ", ")))
  
  year_coverage <- monthly_coverage %>%
    filter(month %in% core_months) %>%
    group_by(year) %>%
    summarise(
      n_months_sampled = sum(sampled),
      n_months_possible = length(core_months),
      pct = round(100 * n_months_sampled / n_months_possible, 1),
      .groups = "drop"
    )
  
  problem_years <- year_coverage %>% filter(pct < 100)
  
  if (nrow(problem_years) > 0) {
    cat("WARNING: Some years have incomplete coverage:\n")
    print(problem_years)
  } else {
    cat("OK: All years have complete coverage for these months.\n")
  }
  
  return(year_coverage)
}

coverage_current <- check_core_months(CORE_MONTHS_CURRENT, "Your current selection")
coverage_datadriven <- check_core_months(CORE_MONTHS_DATADRIVEN, "Data-driven selection")

# --- 7. COMPARE CLUSTERING WITH DIFFERENT MONTH SELECTIONS ---
cat("\n=== SIDE-BY-SIDE: Year coverage comparison ===\n")
comparison <- coverage_current %>%
  rename(current_pct = pct, current_n = n_months_sampled) %>%
  select(year, current_n, current_pct) %>%
  left_join(
    coverage_datadriven %>%
      rename(datadriven_pct = pct, datadriven_n = n_months_sampled) %>%
      select(year, datadriven_n, datadriven_pct),
    by = "year"
  )
print(comparison, n = 25)

