
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
# FORMAT AND PRINT RESULTS (updated)
# =============================================================================

cat("\n=== WILCOXON TESTS: PAIRWISE CLUSTER COMPARISONS ===\n\n")

readable_table <- results %>%
  rowwise() %>%
  mutate(
    `Early C1 vs C2` = format_p(p_1),
    `C2 vs Late C1` = format_p(p_2),
    `Early C1 vs Late C1` = format_p(p_3),
    `C1 comb vs C2` = format_p(p_4)
  ) %>%
  ungroup() %>%
  select(Variable, n_Early, n_C2, n_Late, n_Combined, 
         `Early C1 vs C2`, `C2 vs Late C1`, `Early C1 vs Late C1`, `C1 comb vs C2`)

names(readable_table)[2:5] <- c("n (E-C1)", "n (C2)", "n (L-C1)", "n (C1 comb)")

print(as.data.frame(readable_table), row.names = FALSE)




==========================================================
# LATEX TABLE OUTPUT (updated - removed n_Combined, fits textwidth)
# =============================================================================

latex_table <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Comparison of environmental and biological variables across three clusters: ",
  "Early Cluster~1 (1996-2003), Cluster~2 (2004-2013), and Late Cluster~1 (2014-2016). ",
  "P-values are from two-sided Wilcoxon rank-sum tests on the monthly values. ",
  "Significance levels indicated by asterisks (p-value: ***\\textless0.001, **\\textless0.01, *\\textless0.05), ",
  "non-significant results show the p-value. ",
  "Note that sample sizes for Late Cluster~1 are smaller ($n$ = 21-36), which may limit ",
  "statistical power for comparisons involving this period. ",
  "A pattern consistent with recovery to pre-2004 conditions is indicated when E-C1 vs C2 and ",
  "C2 vs L-C1 show significant differences, while E-C1 vs L-C1 does not.}\n",
  "\\label{tab:ClustCompWilcox}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{l rrr cccc}\n",
  "\\hline\n",
  "& \\multicolumn{3}{c}{Sample size ($n$)} & \\multicolumn{4}{c}{Wilcoxon test} \\\\\n",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-8}\n",
  "Variable & E-C1 & C2 & L-C1 & E-C1 vs C2 & C2 vs L-C1 & E-C1 vs L-C1 & C1 comb vs C2 \\\\\n",
  "\\hline\n"
)

# Add data rows
for (i in 1:nrow(results)) {
  latex_table <- paste0(latex_table,
                        sprintf("%s & %d & %d & %d & %s & %s & %s & %s \\\\\n",
                                results$Variable[i],
                                results$n_Early[i], results$n_C2[i], results$n_Late[i],
                                format_p_latex(results$p_1[i]),
                                format_p_latex(results$p_2[i]),
                                format_p_latex(results$p_3[i]),
                                format_p_latex(results$p_4[i]))
  )
}

# Close table
latex_table <- paste0(latex_table,
                      "\\hline\n",
                      "\\end{tabular}}\n",
                      "\\end{table}\n"
)

cat("\n=== LATEX TABLE ===\n\n")
cat(latex_table)