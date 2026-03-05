library(tidyverse)
require(cowplot)
library(cowplot)
library(ggpubr)

# get env data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

### ADD CLUSTER COMP ###

CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")

# take negative of u10 (vector in u direction, mostly negative) as positive "wind speed":
CARIACO$u10_negative = -CARIACO$u10

# Convert PrimaryProductivity from mg to g C m⁻² d⁻¹
CARIACO$PrimaryProductivity_g <- CARIACO$PrimaryProductivity / 1000

# =============================================================================
# DEFINE 3 CLUSTERS (matching hierarchical clustering from Figure 4)
# =============================================================================

CARIACO <- CARIACO %>%
  mutate(
    year = as.numeric(format(date, "%Y")),
    group = case_when(
      year >= 1996 & year <= 2003 ~ "Early Cluster 1 (1996-2003)",
      year >= 2004 & year <= 2013 ~ "Cluster 2 (2004-2013)",
      year >= 2014 & year <= 2016 ~ "Late Cluster 1 (2014-2016)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group))

# Set factor levels for consistent ordering in plots
CARIACO$group <- factor(CARIACO$group, levels = c("Early Cluster 1 (1996-2003)", 
                                                  "Cluster 2 (2004-2013)", 
                                                  "Late Cluster 1 (2014-2016)"))

# Verify sample sizes
cat("=== SAMPLE SIZES BY CLUSTER ===\n")
CARIACO %>%
  group_by(group) %>%
  summarise(
    n_total = n(),
    years = paste(min(year), max(year), sep = "-"),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# FOLD REDUCTION ANALYSIS (Pre-shift vs Post-shift, kept as original)
# =============================================================================

target_vars <- c("Chlorophyll", 
                 "Abundance_Diatom", 
                 "Abundance_Hapto", 
                 "Abundance_Dino", 
                 "Abundance_Nanoflagellate", 
                 "Abundance_Cyano")

period_comparison <- CARIACO %>%
  mutate(Period_Label = case_when(
    year >= 1996 & year <= 2003 ~ "Period_Early (1996-2003)",
    year >= 2004 & year <= 2013 ~ "Period_Low (2004-2013)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Period_Label))

period_stats <- period_comparison %>%
  group_by(Period_Label) %>%
  summarise(across(all_of(target_vars), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

reduction_analysis <- period_stats %>%
  pivot_longer(cols = -Period_Label, names_to = "Variable", values_to = "Mean_Value") %>%
  pivot_wider(names_from = Period_Label, values_from = Mean_Value) %>%
  rename(
    Baseline = `Period_Early (1996-2003)`,
    Low_Phase = `Period_Low (2004-2013)`
  ) %>%
  mutate(
    Fold_Decrease = Baseline / Low_Phase,
    Pct_Reduction = (Baseline - Low_Phase) / Baseline * 100
  ) %>%
  mutate(across(where(is.numeric), round, 2))

print(reduction_analysis)

# =============================================================================
# PREPARE DATA FOR DENSITY PLOTS
# =============================================================================

ENV_DATA_groups <- CARIACO

envdatmelt <- ENV_DATA_groups %>% 
  select(-date, -time_month, -year) %>% 
  gather(variable, value, -group)

head(envdatmelt)

# =============================================================================
# DENSITY PLOT HELPER FUNCTION
# =============================================================================

# Color palette matching Figure 4
colpallete = c("Early Cluster 1 (1996-2003)" = "#2166AC", 
               "Cluster 2 (2004-2013)" = "#B2182B",
               "Late Cluster 1 (2014-2016)" = "#92C5DE")

make_density_plot <- function(data, var_name, plot_title, x_label, use_log = FALSE, show_y_label = TRUE) {
  
  p <- data %>% 
    filter(variable == var_name) %>%
    ggdensity(x = "value",
              add = "median", rug = TRUE,
              color = "group", fill = "group",
              palette = colpallete) + 
    ggtitle(plot_title) + 
    theme(legend.position = "none") + 
    xlab(x_label)
  
  if (use_log) p <- p + scale_x_log10()
  if (!show_y_label) p <- p + ylab(NULL)
  
  return(p)
}

# =============================================================================
# DEFINE PLOT SPECIFICATIONS
# =============================================================================

plot_specs <- list(
  # Row 1: AMO, MEI v2, Wind speed, Precipitation, Salinity
  A = list(var = "AMO", title = "AMO", xlab = "[AMO]"),
  B = list(var = "MEIv2", title = "MEI v2", xlab = "[MEI v2]", ylab = FALSE),
  C = list(var = "u10_negative", title = "Wind speed", xlab = "[m/s]", ylab = FALSE),
  D = list(var = "tp", title = "Precipitation", xlab = "[m Water/day]", ylab = FALSE),
  E = list(var = "Salinity_bottles", title = "Salinity", xlab = "[PSU]", ylab = FALSE),
  
  # Row 2: SST, 21°C Isotherm, NO3, PO4, SiO4
  F = list(var = "sst_10m", title = "SST", xlab = "[°C]"),
  G = list(var = "Isotherm_21", title = "21°C Isotherm", xlab = "[m]", ylab = FALSE),
  H = list(var = "NO3_merged", title = expression(NO[3]), xlab = "[µM]", ylab = FALSE),
  I = list(var = "PO4_merged", title = expression(PO[4]), xlab = "[µM]", ylab = FALSE),
  J = list(var = "SiO4_merged", title = expression(SiO[4]), xlab = "[µM]", ylab = FALSE),
  
  # Row 3: Primary Productivity, Chlorophyll a, Genus Richness, Shannon Index, Pielou Index
  K = list(var = "PrimaryProductivity_g", title = "Primary Productivity", 
           xlab = expression(paste("[g C ", m^{-2}, " ", d^{-1}, "]")), log = TRUE),
  L = list(var = "Chlorophyll", title = "Chlorophyll a", 
           xlab = expression(paste("[mg ", m^{-2}, "]")), log = TRUE, ylab = FALSE),
  M = list(var = "GenusRichness", title = "Genus Richness", xlab = "[Richness]", ylab = FALSE),
  N = list(var = "Shannon_gen", title = "Shannon Index", xlab = "[Shannon Index]", ylab = FALSE),
  O = list(var = "Pielou_gen", title = "Pielou Index", xlab = "[Pielou Index]", ylab = FALSE)
)

# =============================================================================
# GENERATE ALL PLOTS
# =============================================================================

plots <- lapply(plot_specs, function(spec) {
  make_density_plot(
    data = envdatmelt,
    var_name = spec$var,
    plot_title = spec$title,
    x_label = spec$xlab,
    use_log = ifelse(is.null(spec$log), FALSE, spec$log),
    show_y_label = ifelse(is.null(spec$ylab), TRUE, spec$ylab)
  )
})

# =============================================================================
# CREATE LEGEND
# =============================================================================

legend_plot <- envdatmelt %>% 
  filter(variable == "u10_negative") %>%
  ggdensity(x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + 
  ggtitle("Wind speed") + 
  theme_cowplot(font_size = 20) +
  labs(color = "Periods", fill = "Periods") + 
  theme(legend.position = "top")

legenddd <- get_legend(legend_plot)

# =============================================================================
# COMBINE INTO ROWS
# =============================================================================

prow1 <- plot_grid(plots$A, plots$B, plots$C, plots$D, plots$E, 
                   ncol = 5, labels = c('a', 'b', 'c', 'd', 'e'))

prow2 <- plot_grid(plots$F, plots$G, plots$H, plots$I, plots$J, 
                   ncol = 5, labels = c('f', 'g', 'h', 'i', 'j'))

prow3 <- plot_grid(plots$K, plots$L, plots$M, plots$N, plots$O, 
                   ncol = 5, labels = c('k', 'l', 'm', 'n', 'o'))

# =============================================================================
# FINAL COMBINED PLOT
# =============================================================================

options(repr.plot.width = 12, repr.plot.height = 8)

plottt <- plot_grid(legenddd, prow1, prow2, prow3, 
                    rel_heights = c(0.3, 1, 1, 1), 
                    ncol = 1, align = 'v', axis = 'r')

plottt

#ggsave("plots/exports/Figure5_DensityDistributions_v3.pdf", plottt, width=13, height=7)

# =============================================================================
# PAIRWISE WILCOXON STATISTICS (3 clusters)
# =============================================================================

# Updated variable list including functional groups and Primary Productivity
test_vars = c("AMO", "MEIv2", "u10_negative", "tp", "Salinity_bottles",
              "sst_10m", "Isotherm_21", "NO3_merged", "PO4_merged", "SiO4_merged",
              "PrimaryProductivity_g", "Chlorophyll", 
              "Abundance_Diatom", "Abundance_Hapto", "Abundance_Dino", 
              "Abundance_Cyano", "Abundance_Nanoflagellate",
              "GenusRichness", "Shannon_gen", "Pielou_gen")

test_vars_finalnames = c("AMO", "MEI v.2", "Wind speed", "Precipitation", "Salinity",
                         "SST", "21°C Isotherm", "NO$_3$", "PO$_4$", "SiO$_4$",
                         "Primary Productivity", "Chlorophyll $a$", 
                         "Diatoms", "Haptophytes", "Dinoflagellates", 
                         "Cyanobacteria", "Nanoflagellates",
                         "Genus Richness", "Shannon Index", "Pielou Index")

# Define pairwise comparisons
comparisons <- list(
  c("Early Cluster 1 (1996-2003)", "Cluster 2 (2004-2013)"),
  c("Cluster 2 (2004-2013)", "Late Cluster 1 (2014-2016)"),
  c("Early Cluster 1 (1996-2003)", "Late Cluster 1 (2014-2016)")
)

comparison_names <- c("Early C1 vs C2", "C2 vs Late C1", "Early C1 vs Late C1")

# Run pairwise Wilcoxon tests
results <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  vname <- test_vars_finalnames[i]
  
  # Sample sizes per cluster
  n_early <- sum(!is.na(ENV_DATA_groups[[var]][ENV_DATA_groups$group == "Early Cluster 1 (1996-2003)"]))
  n_c2 <- sum(!is.na(ENV_DATA_groups[[var]][ENV_DATA_groups$group == "Cluster 2 (2004-2013)"]))
  n_late <- sum(!is.na(ENV_DATA_groups[[var]][ENV_DATA_groups$group == "Late Cluster 1 (2014-2016)"]))
  
  row <- data.frame(Variable = vname, n_Early = n_early, n_C2 = n_c2, n_Late = n_late)
  
  for (j in seq_along(comparisons)) {
    comp <- comparisons[[j]]
    
    test_data <- ENV_DATA_groups %>%
      filter(group %in% comp) %>%
      filter(!is.na(.data[[var]]))
    
    if (nrow(test_data) > 0 & length(unique(test_data$group)) == 2) {
      wt <- wilcox.test(as.formula(paste("`", var, "` ~ group", sep = "")), data = test_data)
      row[[paste0("p_", j)]] <- wt$p.value
    } else {
      row[[paste0("p_", j)]] <- NA
    }
  }
  
  results <- rbind(results, row)
}

# =============================================================================
# ADD COMBINED CLUSTER 1 ANALYSIS
# =============================================================================

# Create combined Cluster 1 group
ENV_DATA_groups <- ENV_DATA_groups %>%
  mutate(group_combined = case_when(
    group %in% c("Early Cluster 1 (1996-2003)", "Late Cluster 1 (2014-2016)") ~ "Cluster 1 (combined)",
    group == "Cluster 2 (2004-2013)" ~ "Cluster 2 (2004-2013)",
    TRUE ~ NA_character_
  ))

# Add combined columns to results
for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  
  # Sample size for combined C1
  n_combined <- sum(!is.na(ENV_DATA_groups[[var]][ENV_DATA_groups$group_combined == "Cluster 1 (combined)"]))
  results$n_Combined[i] <- n_combined
  
  # Wilcoxon test: Combined C1 vs C2
  test_data <- ENV_DATA_groups %>%
    filter(!is.na(group_combined)) %>%
    filter(!is.na(.data[[var]]))
  
  if (nrow(test_data) > 0 & length(unique(test_data$group_combined)) == 2) {
    wt <- wilcox.test(as.formula(paste("`", var, "` ~ group_combined", sep = "")), data = test_data)
    results$p_4[i] <- wt$p.value
  } else {
    results$p_4[i] <- NA
  }
}

# =============================================================================
# FDR Correction
# =============================================================================

# Apply FDR correction
all_pvals <- c(results$p_1, results$p_2, results$p_3, results$p_4)
adjusted_pvals <- p.adjust(all_pvals, method = "BH")

# Reshape back
n <- nrow(results)
ç$p_1_adj <- adjusted_pvals[1:n]
results$p_2_adj <- adjusted_pvals[(n+1):(2*n)]
results$p_3_adj <- adjusted_pvals[(2*n+1):(3*n)]
results$p_4_adj <- adjusted_pvals[(3*n+1):(4*n)]


# =============================================================================
# LATEX TABLE WITH ADJUSTED P-VALUES
# =============================================================================

format_p_latex_adj <- function(p_raw, p_adj) {
  # Show significance based on adjusted p-value
  if (is.na(p_adj)) return("-")
  if (p_adj < 0.001) return("***")
  if (p_adj < 0.01) return("**")
  if (p_adj < 0.05) return("*")
  # For non-significant, show the adjusted p-value
  
  return(sprintf("{\\scriptsize %.2f}", p_adj))
}

latex_table <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Comparison of environmental and biological variables across three clusters: ",
  "Early Cluster~1 (1996--2003), Cluster~2 (2004--2013), and Late Cluster~1 (2014--2016). ",
  "P-values are from two-sided Wilcoxon rank-sum tests on monthly values, ",
  "adjusted for multiple comparisons using the Benjamini-Hochberg procedure (80 tests). ",
  "Significance levels: ***$p_{adj} < 0.001$, **$p_{adj} < 0.01$, *$p_{adj} < 0.05$; ",
  "non-significant results show the adjusted p-value. ",
  "Note that sample sizes for Late Cluster~1 are smaller ($n$ = 21--36), which may limit ",
  "statistical power for comparisons involving this period. ",
  "Monthly time series may exhibit temporal autocorrelation, so these tests should be ",
  "interpreted as descriptive comparisons rather than strict hypothesis tests. ",
  "A pattern consistent with recovery to pre-2004 conditions is indicated when E-C1 vs C2 and ",
  "C2 vs L-C1 show significant differences, while E-C1 vs L-C1 does not. ",
  "The combined Cluster~1 (C1) pools Early and Late Cluster~1 data.}\n",
  "\\label{tab:ClustCompWilcox}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l rrr cccc}\n",
  "\\hline\n",
  "& \\multicolumn{3}{c}{Sample size ($n$)} & \\multicolumn{4}{c}{Wilcoxon test ($p_{adj}$)} \\\\\n",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-8}\n",
  "Variable & E-C1 & C2 & L-C1 & E-C1 vs C2 & C2 vs L-C1 & E-C1 vs L-C1 & C1 vs C2 \\\\\n",
  "\\hline\n"
)

# Add data rows with adjusted p-values
for (i in 1:nrow(results)) {
  latex_table <- paste0(latex_table,
                        sprintf("%s & %d & %d & %d & %s & %s & %s & %s \\\\\n",
                                results$Variable[i],
                                results$n_Early[i], results$n_C2[i], results$n_Late[i],
                                format_p_latex_adj(results$p_1[i], results$p_1_adj[i]),
                                format_p_latex_adj(results$p_2[i], results$p_2_adj[i]),
                                format_p_latex_adj(results$p_3[i], results$p_3_adj[i]),
                                format_p_latex_adj(results$p_4[i], results$p_4_adj[i])
                        )
  )
}

latex_table <- paste0(latex_table,
                      "\\hline\n",
                      "\\end{tabular}}\n",
                      "\\end{table}\n"
)

cat(latex_table)




# =============================================================================
# PATTERN SUMMARY
# =============================================================================

cat("\n\n=== PATTERN SUMMARY ===\n")
cat("Variables showing 'recovery' pattern (E-C1 vs C2: sig, C2 vs L-C1: sig, E-C1 vs L-C1: n.s.):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1_adj[i]
  p2 <- results$p_2_adj[i]
  p3 <- results$p_3_adj[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 < 0.05 & p2 < 0.05 & p3 >= 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}

cat("\nVariables showing 'sustained shift' (E-C1 vs C2: sig, C2 vs L-C1: n.s., E-C1 vs L-C1: sig):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1_adj[i]
  p2 <- results$p_2_adj[i]
  p3 <- results$p_3_adj[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 < 0.05 & p2 >= 0.05 & p3 < 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}

cat("\nVariables showing 'partial recovery' (all three comparisons significant):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1_adj[i]
  p2 <- results$p_2_adj[i]
  p3 <- results$p_3_adj[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 < 0.05 & p2 < 0.05 & p3 < 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}

cat("\nVariables showing no regime effect (none significant):\n")

for (i in 1:nrow(results)) {
  p1 <- results$p_1_adj[i]
  p2 <- results$p_2_adj[i]
  p3 <- results$p_3_adj[i]
  
  if (!is.na(p1) & !is.na(p2) & !is.na(p3)) {
    if (p1 >= 0.05 & p2 >= 0.05 & p3 >= 0.05) {
      cat(sprintf("  - %s\n", results$Variable[i]))
    }
  }
}
