# ==============================================================================
# Gradient Forest Analysis of CARIACO Phytoplankton Community
# ==============================================================================
#
# Purpose: Identify environmental drivers of phytoplankton community composition
#          in the CARIACO ocean time series using gradient forest modeling.
#
# Method:  Gradient Forest (Ellis et al. 2012) extends random forest regression
#          to identify nonlinear thresholds in species-environment relationships
#          and quantify predictor importance across the full community.
#
# Outputs:
#   - GF_output_plot1: Combined figure with variable importance and cumulative
#                      importance curves for key predictors
#   - SupplementCumImp_Plot: Cumulative importance for secondary predictors
#   - Actual_Genus.csv: List of genera included in the analysis
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------

# Toggle between lagged and non-lagged environmental predictors
# Set TRUE for analysis with time-lagged climate indices
# Set FALSE for contemporaneous predictor analysis
USE_LAGGED_PREDICTORS <- TRUE

# Species filtering: minimum samples where genus must be present (>0 counts)
# Conservative filter; GF model performs additional internal filtering
MIN_PRESENCE_THRESHOLD <- 10

# Gradient Forest parameters
GF_NTREE <- 1500           # Number of trees per species model
GF_NBIN <- 201             # Bins for cumulative importance curves
GF_CORR_THRESHOLD <- 0.5   # Correlation threshold for predictor selection

# ------------------------------------------------------------------------------
# 2. LOAD LIBRARIES
# ------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(vegan)
library(viridis)
library(cowplot)
library(RColorBrewer)

# gradientForest installation note:
# Requires functional gfortran compiler (install gcc via Homebrew on Mac)
# On R 4.3.3 / macOS Sequoia / M1, may require modified Makeconf
# See: https://stackoverflow.com/questions/76096681/
# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)

# ------------------------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------------------------

# Phytoplankton counts (monthly integrated cell counts, interpolated)
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts), ]

# Standardize genus naming for R compatibility
phyto_counts[phyto_counts$Genus == "Pseudo-nitzschia", ]$Genus <- "Pseudo.nitzschia"

# Environmental data (monthly resolution, pre-computed lags included)
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# convert u10 to Wind Speed to be consistent with manuscript
CARIACO$u10 <- -CARIACO$u10
CARIACO$u10_lag1 <- -CARIACO$u10_lag1
# ------------------------------------------------------------------------------
# 4. DEFINE ENVIRONMENTAL PREDICTOR SETS
# ------------------------------------------------------------------------------

# Lagged predictors: lags determined from preliminary cross-correlation analysis
# Climate indices (AMO, MEI) exhibit multi-year lags in community response
env_factors_lagged <- c(
  "u10_lag1",              # 10m wind speed (1-month lag)
  "NO3_merged",            # Nitrate concentration (contemporaneous)
  "PO4_merged",            # Phosphate concentration
  "SiO4_merged",           # Silicate concentration
  "Salinity_bottles_lag2", # Salinity (2-month lag)
  "sst_10m",               # Sea surface temperature at 10m
  "Isotherm_21",           # 21°C isotherm depth (upwelling proxy)
  "AMO_lag24",             # Atlantic Multidecadal Oscillation (24-month lag)
  "MEIv2_lag48",           # Multivariate ENSO Index v2 (48-month lag)
  "tp_lag1"                # Total precipitation (1-month lag)
)

# Non-lagged predictors for sensitivity analysis
env_factors_nolag <- c(
  "u10",
  "NO3_merged",
  "PO4_merged",
  "SiO4_merged",
  "Salinity_bottles",
  "sst_10m",
  "Isotherm_21",
  "AMO",
  "MEIv2",
  "tp"
)

env_factors <- if (USE_LAGGED_PREDICTORS) env_factors_lagged else env_factors_nolag

# ------------------------------------------------------------------------------
# 5. DATA EXTRACTION FUNCTION
# ------------------------------------------------------------------------------

#' Extract and align phytoplankton counts with environmental data for GF
#'
#' Aggregates species-level counts to genus level, creates month-year key,
#' and performs inner join with environmental data on complete cases only.
#'
#' @param env_factors Character vector of environmental variable names
#' @param phyto_data Phytoplankton count dataframe with Genus, date, counts
#' @param env_data Environmental dataframe with time_month column
#' @return List: env (environmental matrix), genus (abundance matrix)
extract_gf_matrices <- function(env_factors, phyto_data, env_data) {
  
  # Aggregate counts by genus and date
  genus_summary <- phyto_data %>%
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
    group_by(Genus, date) %>%
    summarise(Total = sum(counts), .groups = "drop") %>%
    arrange(desc(date))
  
  # Pivot to wide format (genera as columns, dates as rows)
  genus_wide <- pivot_wider(
    genus_summary,
    names_from = Genus,
    values_from = Total,
    values_fn = sum,
    values_fill = 0
  )
  
  # Create month-year key for joining (format matches env_data)
  genus_wide$time_month <- format(genus_wide$date, format = "%m-%Y")
  genus_wide <- genus_wide %>% select(-date)
  
  # Extract environmental variables
  env_subset <- env_data %>%
    select(time_month, all_of(env_factors))
  
  # Join datasets and remove time key
  
  combined <- list(env_subset, genus_wide) %>%
    reduce(full_join, by = "time_month") %>%
    select(-time_month)
  
  # Retain only complete cases (listwise deletion)
  complete_data <- combined[complete.cases(combined), ]
  
  # Separate matrices
  env_matrix <- complete_data %>% select(all_of(env_factors))
  genus_matrix <- complete_data %>% select(-all_of(env_factors))
  
  return(list(env = env_matrix, genus = genus_matrix))
}

# ------------------------------------------------------------------------------
# 6. PREPARE DATA FOR GRADIENT FOREST
# ------------------------------------------------------------------------------

# Extract aligned matrices
gf_data <- extract_gf_matrices(env_factors, phyto_counts, CARIACO)
env_matrix <- gf_data$env
genus_matrix <- gf_data$genus

# Filter rare genera (present in fewer than threshold samples)
genus_filtered <- genus_matrix %>%
  select_if(colSums(. > 0) > MIN_PRESENCE_THRESHOLD)

# Clean column names for R compatibility
colnames(genus_filtered) <- make.names(colnames(genus_filtered))

# Log-transform abundance data
# log1p(x) = log(1 + x) handles zeros appropriately
genus_log <- log1p(genus_filtered)

# Calculate optimal tree depth following Pitcher et al. (2012)
# maxLevel controls tree complexity to avoid overfitting
n_sites <- nrow(genus_filtered)
n_species <- ncol(genus_filtered)
max_level <- floor(log2(n_sites * 0.368 / 2))

cat("=== Gradient Forest Input Summary ===\n")
cat("  Time points (sites):", n_sites, "\n")
cat("  Genera (species):", n_species, "\n")
cat("  Max tree level:", max_level, "\n")
cat("  Predictor set:", ifelse(USE_LAGGED_PREDICTORS, "LAGGED", "NO-LAG"), "\n\n")

# ------------------------------------------------------------------------------
# 7. FIT GRADIENT FOREST MODEL
# ------------------------------------------------------------------------------

cat("Fitting Gradient Forest model...\n")

gf_model <- gradientForest(
  cbind(env_matrix, genus_log),
  predictor.vars = colnames(env_matrix),
  response.vars = colnames(genus_filtered),
  ntree = GF_NTREE,
  transform = NULL,       # Data already log-transformed
  compact = TRUE,
  trace = TRUE,
  nbin = GF_NBIN,
  maxLevel = max_level,
  corr.threshold = GF_CORR_THRESHOLD
)

print(gf_model)

# ------------------------------------------------------------------------------
# 8. EXTRACT IMPORTANCE METRICS
# ------------------------------------------------------------------------------

# Extract top 4 predictors by weighted importance for main figure
main_predictors <- names(importance(gf_model))[1:4]
cat("Main predictors for figure:", paste(main_predictors, collapse = ", "), "\n")

# Weighted R² importance for predictors
imp_weighted <- importance(gf_model, "Weighted")
imp_order <- order(imp_weighted)

importance_df <- data.frame(
  variable = factor(names(imp_weighted), levels = names(imp_weighted)[imp_order]),
  importance = imp_weighted,
  row.names = NULL
)

# Species-level R² importance
species_importance <- importance(gf_model, type = "Species")
species_order <- order(species_importance)

species_imp_df <- data.frame(
  Genus = factor(names(species_importance)[species_order],
                 levels = names(species_importance)[species_order]),
  importance = species_importance[species_order],
  row.names = NULL
)

# ------------------------------------------------------------------------------
# 9. EXTRACT CUMULATIVE IMPORTANCE CURVES
# ------------------------------------------------------------------------------

imp_vars <- names(imp_weighted)

# --- Species-level cumulative importance ---
species_cumimp_list <- list()

for (var in imp_vars) {
  CU <- cumimp(gf_model, var, "Species")
  
  for (sp in names(CU)) {
    # Subsample points for plotting efficiency
    n_pts <- length(CU[[sp]]$x)
    isub <- seq(1, n_pts, length.out = min(500, n_pts))
    
    species_cumimp_list[[paste(var, sp, sep = "_")]] <- data.frame(
      Predictor = var,
      Species = sp,
      x = CU[[sp]]$x[isub],
      y = CU[[sp]]$y[isub]
    )
  }
}

species_cumimp_df <- rbindlist(species_cumimp_list)
species_cumimp_df$Predictor <- factor(species_cumimp_df$Predictor, levels = names(imp_weighted))

# --- Overall cumulative importance ---
overall_cumimp_list <- list()

for (var in imp_vars) {
  CU <- cumimp(gf_model, var)
  n_pts <- length(CU$x)
  isub <- seq(1, n_pts, length.out = min(500, n_pts))
  
  overall_cumimp_list[[var]] <- data.frame(
    Predictor = var,
    x = CU$x[isub],
    y = CU$y[isub]
  )
}

overall_cumimp_df <- rbindlist(overall_cumimp_list)
overall_cumimp_df$Predictor <- factor(overall_cumimp_df$Predictor, levels = names(imp_weighted))

# ------------------------------------------------------------------------------
# 10. FUNCTIONAL GROUP CUMULATIVE IMPORTANCE
# ------------------------------------------------------------------------------

# Get genera actually modeled (extracted from species cumimp, mirrors original approach)
modeled_genera <- unique(species_cumimp_df$Species)

# Match genera to functional groups
# Apply make.names() to phyto_counts$Genus to ensure names match model output
genus_funcgroup <- phyto_counts %>%
  mutate(Genus = make.names(Genus)) %>%
  group_by(Genus) %>%
  summarise(FuncGroup = first(FuncGroup), .groups = "drop")

genus_matched <- genus_funcgroup %>%
  filter(Genus %in% modeled_genera)

# Report taxonomic coverage
coverage <- sum(genus_matrix[, genus_matched$Genus]) / sum(genus_matrix)
cat("\nProportion of total counts in model:", round(coverage, 3), "\n")

# Define colors for Functional Groups - hardcoded for consistency across scripts
funcgroup_colors <- c(
  "Dinoflagellata" = "#E41A1C",
  "Bacillariophyceae" = "#377EB8",
  "Haptophyta" = "#4DAF4A",
  "Cyanobacteria" = "#984EA3",
  "Others" = "#FF7F00"
)

genus_matched$color <- funcgroup_colors[genus_matched$FuncGroup]
species_colors <- setNames(genus_matched$color, genus_matched$Genus)
# --- Helper functions for group-level cumulative importance ---
# Adapted from gradientForest source code

inverse_density <- function(dens) {
  dens$y <- 1 / dens$y
  dens
}

normalize_density <- function(f) {
  integral <- try(integrate(approxfun(f, rule = 2),
                            lower = min(f$x), upper = max(f$x))$value,
                  silent = TRUE)
  if (inherits(integral, "try-error")) {
    integral <- sum(f$y) * diff(f$x)[1]
  }
  f$y <- f$y / integral * diff(range(f$x))
  f
}

aggregate_sum <- function(x, by, sort.it = FALSE) {
  if (!is.data.frame(by)) by <- data.frame(by)
  if (sort.it) {
    ord <- do.call("order", unname(by))
    x <- x[ord]
    by <- by[ord, , drop = FALSE]
  }
  logical_diff <- function(group) group[-1] != group[-length(group)]
  change <- logical_diff(by[[1]])
  for (i in seq_along(by)[-1]) {
    change <- change | logical_diff(by[[i]])
  }
  by <- by[c(TRUE, change), , drop = FALSE]
  by$x <- diff(c(0, cumsum(x)[c(change, TRUE)]))
  by
}

#' Calculate cumulative importance for a taxonomic group
#'
#' @param importance_df Subset of gf_model$res for one predictor
#' @param rsq Mean R² across species in the group
#' @param predictor Predictor variable name
#' @param gf_model Gradient forest model object
#' @param standardize Standardize by inverse predictor density
#' @param standardize_after Apply standardization after aggregation
#' @return List with x (predictor values) and y (cumulative importance)
get_group_cumimp <- function(importance_df, rsq, predictor, gf_model,
                             standardize = TRUE, standardize_after = FALSE) {
  if (nrow(importance_df) == 0) {
    return(list(x = 0, y = 0))
  }
  
  agg <- with(importance_df, aggregate_sum(improve.norm, list(split), sort.it = TRUE))
  cum_split <- agg[, 1]
  height <- agg[, 2]
  
  if (standardize & standardize_after) {
    dinv <- normalize_density(inverse_density(gf_model$dens[[predictor]]))
  } else {
    dinv <- inverse_density(gf_model$dens[[predictor]])
  }
  
  dinv_vals <- approx(dinv$x, dinv$y, cum_split, rule = 2)$y
  
  if (any(bad <- is.na(height))) {
    cum_split <- cum_split[!bad]
    height <- height[!bad]
    dinv_vals <- dinv_vals[!bad]
  }
  
  if (standardize & !standardize_after) height <- height * dinv_vals
  height <- height / sum(height) * rsq
  if (standardize & standardize_after) height <- height * dinv_vals
  
  list(x = cum_split, y = cumsum(height))
}

# --- Calculate functional group cumulative importance ---
funcgroup_cumimp_list <- list()

for (var in imp_vars) {
  importance_subset <- gf_model$res[gf_model$res$var == var, ]
  
  for (fg in unique(genus_matched$FuncGroup)) {
    genera_in_group <- genus_matched %>%
      filter(FuncGroup == fg) %>%
      pull(Genus)
    
    fg_cumimp <- get_group_cumimp(
      subset(importance_subset, spec %in% genera_in_group),
      mean(gf_model$imp.rsq[var, genera_in_group]),
      var,
      gf_model
    )
    
    n_pts <- length(fg_cumimp$x)
    isub <- seq(1, n_pts, length.out = min(500, n_pts))
    
    funcgroup_cumimp_list[[paste(var, fg, sep = "_")]] <- data.frame(
      Predictor = var,
      FuncGroup = fg,
      x = fg_cumimp$x[isub],
      y = fg_cumimp$y[isub]
    )
  }
}

funcgroup_cumimp_df <- rbindlist(funcgroup_cumimp_list)
funcgroup_cumimp_df$Predictor <- factor(funcgroup_cumimp_df$Predictor,
                                        levels = levels(overall_cumimp_df$Predictor))

# ------------------------------------------------------------------------------
# 11. CREATE PLOTS
# ------------------------------------------------------------------------------

importance_title <- expression(paste(R^2, " weighted importance"))

# --- Predictor importance barplot ---
weighted_imp_plot <- ggplot(importance_df, aes(x = variable, y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  #scale_x_discrete(limits = rev) +
  labs(x = "Predictor Variable", y = importance_title) +
  theme_cowplot()

# --- Species importance dotplot ---
species_imp_plot <- ggplot(species_imp_df, aes(x = Genus, y = importance, color = Genus)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(y = importance_title) +
  scale_color_manual(values = species_colors) +
  guides(color = "none") +
  theme_cowplot(font_size = 10)

# --- Functional group legend ---
legend_data <- data.frame(
  FuncGroup = names(funcgroup_colors),
  y = seq_along(funcgroup_colors)
)

legend_plot <- ggplot(legend_data, aes(x = FuncGroup, y = y, color = FuncGroup)) +
  geom_point() +
  scale_color_manual(values = funcgroup_colors) +
  theme_cowplot()

funcgroup_legend <- get_legend(legend_plot)

# --- Main cumulative importance plot (key predictors) ---
main_cumimp_plot <- ggplot(funcgroup_cumimp_df %>% filter(Predictor %in% main_predictors)) +
  geom_line(aes(x = x, y = y, color = FuncGroup), linewidth = 1.1) +
  geom_line(
    data = overall_cumimp_df %>% filter(Predictor %in% main_predictors),
    aes(x = x, y = y),
    linewidth = 1.5,
    linetype = "22"
  ) +
  facet_wrap(~Predictor, scales = "free_x") +
  labs(x = "Predictor Value", y = "Cumulative Importance") +
  scale_color_manual(values = funcgroup_colors) +
  guides(color = "none") +
  theme_cowplot(font_size = 20)

# --- Supplemental cumulative importance plot (secondary predictors) ---
SupplementCumImp_Plot <- ggplot(funcgroup_cumimp_df %>% filter(!Predictor %in% main_predictors)) +
  geom_line(aes(x = x, y = y, color = FuncGroup), linewidth = 1.1,alpha=0.7) +
  geom_line(
    data = overall_cumimp_df %>% filter(!Predictor %in% main_predictors),
    aes(x = x, y = y),
    linewidth = 1.5,
    linetype = "22"
  ) +
  facet_wrap(~Predictor, scales = "free_x") +
  labs(x = "Predictor Value", y = "Cumulative Importance") +
  scale_color_manual(values = funcgroup_colors) +
  guides(color = "none") +
  theme_cowplot(font_size = 20)

# ------------------------------------------------------------------------------
# 12. ASSEMBLE COMBINED FIGURE
# ------------------------------------------------------------------------------

right_column <- plot_grid(
  weighted_imp_plot,
  species_imp_plot,
  funcgroup_legend,
  ncol = 3,
  rel_widths = c(1, 1, 0.1)
)

GF_output_plot1 <- plot_grid(
  right_column,
  main_cumimp_plot,
  ncol = 1,
  rel_heights = c(1.8, 2)
)

# Display plots
print(GF_output_plot1)
print(SupplementCumImp_Plot)

# ------------------------------------------------------------------------------
# 13. EXPORT OUTPUTS
# ------------------------------------------------------------------------------

# Determine output suffix based on analysis type
output_suffix <- ifelse(USE_LAGGED_PREDICTORS, "LAGGED", "NOLAG")

# save figures:
ggsave(
  paste0("plots/exports/Figure6_NEW_GF_output_", output_suffix, ".pdf"),
  GF_output_plot1,
  width = 12, height = 16
)
ggsave(
  paste0("plots/exports/Figure6_NEW_GF_supplement_", output_suffix, ".pdf"),
  SupplementCumImp_Plot,
  width = 12, height = 8
)








# Export genus list for manuscript
input_genus <- names(mo5_specGF)
#GenusInput <- GenusFuncGroup_match[GenusFuncGroup_match$Genus %in% input_genus,]
#GenusInput %>% arrange(FuncGroup) %>% select(FuncGroup, Genus) %>% kbl("latex", booktabs=T)

ds_genus_overview <- phyto_counts %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species")%>%
  group_by(Genus, FuncGroup) %>% summarize(IDS = unique(AphiaID))

Actual_Genus = ds_genus_overview[ds_genus_overview$Genus %in% input_genus,]

Actual_Genus$name <- NA
Actual_Genus$auth <- NA
Actual_Genus$rank <- NA

for (i in 1:nrow(Actual_Genus)){
  #print(i)
  id = Actual_Genus$IDS[i]
  
  skip_to_next <- FALSE
  #print(id)
  tryCatch({
    record = wm_record(id=id)
    #print(record)
    auth = record$authority
    name = record$scientificname
    rank = record$rank
    
    Actual_Genus$auth[i] <- auth
    Actual_Genus$name[i] <- name
    Actual_Genus$rank[i] <- rank
    
  }, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
}

#write.csv(Actual_Genus, "Actual_Genus_FULL.csv")