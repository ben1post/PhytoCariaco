# ============================================
# FULL SENSITIVITY ANALYSIS TABLE
# ============================================

library(tidyverse)

# --- 1. LOAD AND PREP DATA ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    year = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m")),
    u10_negative = -u10,
    period = case_when(
      year >= 1995 & year <= 2003 ~ "Early_C1",
      year >= 2004 & year <= 2013 ~ "C2",
      year >= 2014 & year <= 2017 ~ "Late_C1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period))

# --- 2. DEFINE VARIABLES (untransformed) ---
vars <- c("MEIv2", "AMO", "u10_negative", "tp", "Isotherm_21", "Salinity_bottles", "sst_10m",
          "NO3_merged", "PO4_merged", "SiO4_merged", "Chlorophyll",
          "Abundance_Diatom", "Abundance_Hapto", "Abundance_Dino", 
          "Abundance_Cyano", "Abundance_Nanoflagellate",
          "GenusRichness", "Shannon_gen", "Pielou_gen")

renamed_vars <- c("MEI v.2", "AMO", "Wind speed", "Precipitation", "Isotherm Depth", 
                  "Salinity", "SST", "NO$_3$", "PO$_4$", "SiO$_4$", "Chlorophyll $a$",
                  "Diatoms", "Haptophytes", "Dinoflagellates", "Cyanobacteria", 
                  "Nanoflagellates", "Richness", "Shannon", "Pielou")

# --- 3. VERIFY SAMPLE SIZES ---
cat("=== SAMPLE SIZES BY PERIOD ===\n")
CARIACO %>%
  group_by(period) %>%
  summarise(
    n_total = n(),
    years = paste(min(year), max(year), sep = "-"),
    .groups = "drop"
  ) %>%
  print()

# --- 4. RUN TESTS ---
comparisons <- list(
  c("Early_C1", "C2"),
  c("C2", "Late_C1"),
  c("Early_C1", "Late_C1")
)

results <- data.frame()

for (i in seq_along(vars)) {
  var <- vars[i]
  vname <- renamed_vars[i]
  
  # Sample sizes per period
  n_early <- sum(!is.na(CARIACO[[var]][CARIACO$period == "Early_C1"]))
  n_c2 <- sum(!is.na(CARIACO[[var]][CARIACO$period == "C2"]))
  n_late <- sum(!is.na(CARIACO[[var]][CARIACO$period == "Late_C1"]))
  
  row <- data.frame(Variable = vname, n_Early = n_early, n_C2 = n_c2, n_Late = n_late)
  
  for (j in seq_along(comparisons)) {
    comp <- comparisons[[j]]
    
    test_data <- CARIACO %>%
      filter(period %in% comp) %>%
      filter(!is.na(.data[[var]]))
    
    if (nrow(test_data) > 0 & length(unique(test_data$period)) == 2) {
      wt <- wilcox.test(as.formula(paste(var, "~ period")), data = test_data)
      row[[paste0("p_", j)]] <- wt$p.value
    } else {
      row[[paste0("p_", j)]] <- NA
    }
  }
  
  results <- rbind(results, row)
}

# --- 5. PRINT READABLE TABLE ---
cat("\n=== WILCOXON TESTS: PERIOD COMPARISONS ===\n\n")

format_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return(sprintf("%.2f", p))
}

readable_table <- results %>%
  rowwise() %>%
  mutate(
    `Early C1 vs C2` = format_p(p_1),
    `C2 vs Late C1` = format_p(p_2),
    `Early C1 vs Late C1` = format_p(p_3)
  ) %>%
  ungroup() %>%
  select(Variable, n_Early, n_C2, n_Late, `Early C1 vs C2`, `C2 vs Late C1`, `Early C1 vs Late C1`)

names(readable_table)[2:4] <- c("n (E-C1)", "n (C2)", "n (L-C1)")

print(as.data.frame(readable_table), row.names = FALSE)

# --- 6. LATEX OUTPUT ---
cat("\n\n=== LATEX TABLE ===\n\n")

# Helper function for LaTeX formatting
format_p_latex <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("$^{***}$")
  if (p < 0.01) return("$^{**}$")
  if (p < 0.05) return("$^{*}$")
  return(sprintf("{\\scriptsize %.2f}", p))
}

cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Comparison of environmental and biological variables across three periods: Early Cluster~1 (1995--2004), Cluster~2 (2005--2013), and Late Cluster~1 (2014--2017). P-values are from two-sided Wilcoxon rank-sum tests. Significance levels indicated by asterisks: $^{***}p < 0.001$, $^{**}p < 0.01$, $^{*}p < 0.05$; non-significant results show the p-value. A pattern consistent with partial recovery to pre-2005 conditions is indicated when Early~C1 vs C2 and C2 vs Late~C1 show significant differences, while Early~C1 vs Late~C1 does not.}\n")
cat("\\label{tab:period_comparison}\n")
cat("\\small\n")
cat("\\begin{tabular}{l rrr ccc}\n")
cat("\\hline\n")
cat("& \\multicolumn{3}{c}{Sample size ($n$)} & \\multicolumn{3}{c}{Wilcoxon test} \\\\\n")
cat("\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n")
cat("Variable & E-C1 & C2 & L-C1 & E-C1 vs C2 & C2 vs L-C1 & E-C1 vs L-C1 \\\\\n")
cat("\\hline\n")

for (i in 1:nrow(results)) {
  cat(sprintf("%s & %d & %d & %d & %s & %s & %s \\\\\n",
              results$Variable[i],
              results$n_Early[i], results$n_C2[i], results$n_Late[i],
              format_p_latex(results$p_1[i]),
              format_p_latex(results$p_2[i]),
              format_p_latex(results$p_3[i])))
}

cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")

# --- 7. SUMMARY OF PATTERNS ---
cat("\n\n=== PATTERN SUMMARY ===\n")
cat("Variables showing 'recovery' pattern (E-C1 vs C2: sig, C2 vs L-C1: sig, E-C1 vs L-C1: n.s.):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1[i]
  p2 <- results$p_2[i]
  p3 <- results$p_3[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 < 0.05 & p2 < 0.05 & p3 >= 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}

cat("\nVariables showing 'partial recovery' (E-C1 vs C2: sig, C2 vs L-C1: sig, E-C1 vs L-C1: sig):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1[i]
  p2 <- results$p_2[i]
  p3 <- results$p_3[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 < 0.05 & p2 < 0.05 & p3 < 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}

cat("\nVariables showing no regime effect (none significant):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1[i]
  p2 <- results$p_2[i]
  p3 <- results$p_3[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 >= 0.05 & p2 >= 0.05 & p3 >= 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}