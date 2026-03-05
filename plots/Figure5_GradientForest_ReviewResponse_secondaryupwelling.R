##### Prelim Runs #####


library(tidyverse)
library(vegan)
library(viridis)
library(cowplot)
library(gradientForest)

# Read data
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts),]
phyto_counts[phyto_counts$Genus=="Pseudo-nitzschia",]$Genus <- "Pseudo.nitzschia"

CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# Season function
get_season_3cat <- function(months) {
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  season_lookup[months]
}

# Modified extraction function with season filter
extractMatrixFix_seasonal <- function(env_factors, season_filter = NULL) {
  ds_genus <- phyto_counts %>% 
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts), .groups = "drop") %>%
    arrange(desc(date))
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  Mesh_genus$time_month <- format(Mesh_genus$date, format = "%m-%Y")
  new_Mesh_genus <- Mesh_genus %>% select(-date)
  
  Env_dat <- CARIACO %>% select("time_month", all_of(env_factors))
  
  CARIACO_dat_joined <- list(Env_dat, new_Mesh_genus) %>% 
    reduce(full_join, by = c("time_month"))
  
  if (!is.null(season_filter)) {
    CARIACO_dat_joined <- CARIACO_dat_joined %>%
      mutate(month = substr(time_month, 1, 2),
             season = get_season_3cat(month)) %>%
      filter(season == season_filter) %>%
      select(-month, -season)
  }
  
  CARIACO_dat_joined <- CARIACO_dat_joined %>% select(-time_month)
  Comp_Matrix <- CARIACO_dat_joined[complete.cases(CARIACO_dat_joined),]
  
  Matrix_Env <- Comp_Matrix %>% select(all_of(env_factors))
  Matrix_Genus <- Comp_Matrix %>% select(-all_of(env_factors))
  
  return(list(Matrix_Env, Matrix_Genus))
}

# Run GF model - matches your original runGFmodel()
runGFmodel_seasonal <- function(env_factors, season_filter = NULL) {
  data <- extractMatrixFix_seasonal(env_factors, season_filter)
  
  envGF <- data[[1]]
  specGF <- data[[2]]
  
  mo5_specGF <- specGF %>% select_if(colSums(. > 0) > 10)
  colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))
  
  nSites <- dim(mo5_specGF)[1]
  lev <- floor(log2(nSites * 0.368/2))
  
  log_specGF <- log1p(mo5_specGF * 100)
  
  gf <- gradientForest(cbind(envGF, log_specGF),
                       predictor.vars = colnames(envGF),
                       response.vars = colnames(mo5_specGF),
                       ntree = 1500, transform = NULL, compact = TRUE,
                       trace = TRUE, nbin = 201, maxLevel = lev,
                       corr.threshold = 0.5)
  return(gf)
}

# Define all lag vectors (exactly as in your original)
GF_env_factors_NO3 <- c("NO3_merged", "NO3_merged_lag1", "NO3_merged_lag2", "NO3_merged_lag3")
GF_env_factors_PO4 <- c("PO4_merged", "PO4_merged_lag1", "PO4_merged_lag2", "PO4_merged_lag3")
GF_env_factors_SiO4 <- c("SiO4_merged", "SiO4_merged_lag1", "SiO4_merged_lag2", "SiO4_merged_lag3")
GF_env_factors_SST <- c("sst_10m", "sst_10m_lag1", "sst_10m_lag2", "sst_10m_lag3")
GF_env_factors_Isotherm_21 <- c("Isotherm_21", "Isotherm_21_lag1", "Isotherm_21_lag2", "Isotherm_21_lag3")
GF_env_factors_Salinity <- c("Salinity_bottles", "Salinity_bottles_lag1", "Salinity_bottles_lag2", "Salinity_bottles_lag3")
GF_env_factors_u10 <- c("u10", "u10_lag1", "u10_lag2", "u10_lag3")
GF_env_factors_tp <- c("tp", "tp_lag1", "tp_lag2", "tp_lag3")
GF_env_factors_e <- c("e", "e_lag1", "e_lag2", "e_lag3")
GF_env_factors_AMO <- c("AMO", "AMO_lag1", "AMO_lag2", "AMO_lag3", "AMO_lag4", "AMO_lag5", "AMO_lag6", "AMO_lag12")
GF_env_factors_MEIv2 <- c("MEIv2", "MEIv2_lag1", "MEIv2_lag2", "MEIv2_lag3", "MEIv2_lag4", "MEIv2_lag5", "MEIv2_lag6", "MEIv2_lag12")


# Run lag analysis for each season with error handling
seasons <- c("Upwelling", "Secondary Upwelling", "Rainy")

for (s in seasons) {
  cat("\n\n========================================\n")
  cat("Running lag analysis for:", s, "\n")
  cat("========================================\n")
  
  # Safe wrapper for model runs
  safe_run <- function(env_factors, label) {
    tryCatch({
      gf <- runGFmodel_seasonal(env_factors, s)
      return(as.numeric(importance(gf, sort = FALSE)))
    }, error = function(e) {
      cat("ERROR for", label, ":", e$message, "\n")
      return(rep(NA, length(env_factors)))
    })
  }
  
  NO3_lags <- safe_run(GF_env_factors_NO3, "NO3")
  PO4_lags <- safe_run(GF_env_factors_PO4, "PO4")
  SiO4_lags <- safe_run(GF_env_factors_SiO4, "SiO4")
  SST_lags <- safe_run(GF_env_factors_SST, "SST")
  Isotherm_21_lags <- safe_run(GF_env_factors_Isotherm_21, "Isotherm_21")
  Salinity_lags <- safe_run(GF_env_factors_Salinity, "Salinity")
  u10_lags <- safe_run(GF_env_factors_u10, "u10")
  tp_lags <- safe_run(GF_env_factors_tp, "tp")
  e_lags <- safe_run(GF_env_factors_e, "e")
  AMO_lags <- safe_run(GF_env_factors_AMO, "AMO")
  MEIv2_lags <- safe_run(GF_env_factors_MEIv2, "MEIv2")
  
  # Build tables
  sq <- seq(4)
  sq_climate <- seq(8)
  
  variable_selection <- data.frame(
    "lag" = c(0:3),
    "NO3" = NO3_lags[sq], 
    "PO4" = PO4_lags[sq], 
    "SiO4" = SiO4_lags[sq], 
    "Temperature" = SST_lags[sq],               
    "Isotherm_21" = Isotherm_21_lags[sq], 
    "Salinity" = Salinity_lags[sq],
    "Wind_speed" = u10_lags[sq],
    "Precipitation" = tp_lags[sq],
    "Evaporation" = e_lags[sq]
  )
  
  climate_selection <- data.frame(
    "lag" = c(0:6, 12),
    "AMO" = AMO_lags[sq_climate],
    "MEIv2" = MEIv2_lags[sq_climate]
  )
  
  cat("\n--- In-situ variables (lag 0-3) ---\n")
  print(round(variable_selection, 4))
  
  cat("\n--- Climate variables (lag 0-6, 12) ---\n")
  print(round(climate_selection, 4))
  
  # Save results
  season_label <- gsub(" ", "_", s)
  saveRDS(list(insitu = variable_selection, climate = climate_selection), 
          paste0("GF_lag_results_", season_label, ".RDS"))
}



##### Actual Model Runs #####



library(tidyverse)
library(data.table)
library(vegan)
library(viridis)
library(cowplot)
library(gradientForest)

# Read Phytoplankton Data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts),]
phyto_counts[phyto_counts$Genus=="Pseudo-nitzschia",]$Genus <- "Pseudo.nitzschia"

# Get environmental data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# Define season function
get_season_3cat <- function(months) {
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  season_lookup[months]
}

# Modified extraction function with season filter
extractMatrixFix_seasonal <- function(env_factors, season_filter = NULL) {
  ds_genus <- phyto_counts %>% 
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts)) %>%
    arrange(desc(date))
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  Mesh_genus$time_month <- format(Mesh_genus$date, format = "%m-%Y")
  new_Mesh_genus <- Mesh_genus %>% select(-date)
  
  Env_dat <- CARIACO %>% 
    select("time_month", all_of(env_factors))
  
  CARIACO_dat_joined <- list(Env_dat, new_Mesh_genus) %>% 
    reduce(full_join, by = c("time_month"))
  
  # Apply season filter if specified
  if (!is.null(season_filter)) {
    CARIACO_dat_joined <- CARIACO_dat_joined %>%
      mutate(month = substr(time_month, 1, 2),
             season = get_season_3cat(month)) %>%
      filter(season == season_filter) %>%
      select(-month, -season)
  }
  
  CARIACO_dat_joined <- CARIACO_dat_joined %>% select(-time_month)
  
  Comp_Matrix <- CARIACO_dat_joined[complete.cases(CARIACO_dat_joined),]
  
  Matrix_Env <- Comp_Matrix %>% select(all_of(env_factors))
  Matrix_Genus <- Comp_Matrix %>% select(-all_of(env_factors))
  
  return(list(Matrix_Env, Matrix_Genus))
}

# Non-lagged environmental factors
GF_env_factors_nolag <- c("u10", "NO3_merged", "PO4_merged", "SiO4_merged", 
                          "tp", "e", "Salinity_bottles", "sst_10m", 
                          "Isotherm_21", "AMO", "MEIv2")

# Extract data for SECONDARY UPWELLING only
GF_inputs_secondary <- extractMatrixFix_seasonal(GF_env_factors_nolag, season_filter = "Secondary Upwelling")

envGF <- GF_inputs_secondary[[1]]
specGF <- GF_inputs_secondary[[2]]

# Check sample size
cat("Sample size (Secondary Upwelling):", nrow(envGF), "\n")

# Filter genera present in >10 samples
mo5_specGF <- specGF %>% select_if(colSums(. > 0) > 10)
colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))

cat("Number of genera after filtering:", ncol(mo5_specGF), "\n")

nSites <- dim(mo5_specGF)[1]
nSpecs <- dim(mo5_specGF)[2]
lev <- floor(log2(nSites * 0.368/2))

cat("nSites:", nSites, "\n")
cat("nSpecs:", nSpecs, "\n")
cat("maxLevel:", lev, "\n")

# Log transform species counts
log_specGF <- log1p(mo5_specGF * 100)

# Run Gradient Forest
gf_secondary <- gradientForest(cbind(envGF, log_specGF),
                               predictor.vars = colnames(envGF),
                               response.vars = colnames(mo5_specGF),
                               ntree = 1500,
                               transform = NULL,
                               compact = TRUE,
                               trace = TRUE,
                               nbin = 201,
                               maxLevel = lev,
                               corr.threshold = 0.5)

gf_secondary

# Plot overall importance
plot(gf_secondary, plot.type = "O")

# Top 4 important variables
most_important <- names(importance(gf_secondary))[1:4]
most_important

# Cumulative importance plot
plot(gf_secondary, plot.type = "C", imp.vars = most_important,
     show.species = FALSE, common.scale = TRUE, cex.axis = 0.6, cex.lab = 0.7,
     line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))





###### FINAL RUNS SEASONAL COMP #####

library(tidyverse)
library(data.table)
library(vegan)
library(viridis)
library(cowplot)
library(gradientForest)

# Read data
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts),]
phyto_counts[phyto_counts$Genus=="Pseudo-nitzschia",]$Genus <- "Pseudo.nitzschia"

CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# Season function
get_season_3cat <- function(months) {
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  season_lookup[months]
}

# Extraction function with season filter
extractMatrixFix_seasonal <- function(env_factors, season_filter = NULL) {
  ds_genus <- phyto_counts %>% 
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts), .groups = "drop") %>%
    arrange(desc(date))
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  Mesh_genus$time_month <- format(Mesh_genus$date, format = "%m-%Y")
  new_Mesh_genus <- Mesh_genus %>% select(-date)
  
  Env_dat <- CARIACO %>% select("time_month", all_of(env_factors))
  
  CARIACO_dat_joined <- list(Env_dat, new_Mesh_genus) %>% 
    reduce(full_join, by = c("time_month"))
  
  if (!is.null(season_filter)) {
    CARIACO_dat_joined <- CARIACO_dat_joined %>%
      mutate(month = substr(time_month, 1, 2),
             season = get_season_3cat(month)) %>%
      filter(season == season_filter) %>%
      select(-month, -season)
  }
  
  CARIACO_dat_joined <- CARIACO_dat_joined %>% select(-time_month)
  Comp_Matrix <- CARIACO_dat_joined[complete.cases(CARIACO_dat_joined),]
  
  Matrix_Env <- Comp_Matrix %>% select(all_of(env_factors))
  Matrix_Genus <- Comp_Matrix %>% select(-all_of(env_factors))
  
  return(list(Matrix_Env, Matrix_Genus))
}

# Define optimal lag variable sets based on lag analysis results
# Upwelling: optimal lags from results
GF_env_Upwelling <- c(
  "NO3_merged",           # lag 0
  
  "PO4_merged",           # lag 0
  "SiO4_merged",          # lag 0
  "sst_10m",              # lag 0
  "Isotherm_21",          # lag 0
  "Salinity_bottles",     # lag 0
  "u10",                  # lag 0
  "tp_lag1",              # lag 1
  "AMO_lag3",             # lag 4
  "MEIv2_lag4"            # lag 4
)

# Secondary Upwelling: optimal lags (excluding Isotherm_21 and e which failed)
GF_env_SecondaryUpwelling <- c(
  "NO3_merged",           # lag 0
  "PO4_merged",           # lag 0
  "SiO4_merged_lag1",     # lag 1
  "sst_10m",         # lag 2
  "Isotherm_21",          # lag 0
  "Salinity_bottles_lag2",# lag 2
  "u10",                  # lag 0
  "tp_lag2",              # lag 2
  "AMO_lag2",             # lag 2
  "MEIv2_lag5"            # lag 5
)

# Rainy: optimal lags
GF_env_Rainy <- c(
  "NO3_merged",           # lag 0
  "PO4_merged",           # lag 0
  "SiO4_merged_lag2",     # lag 2
  "sst_10m",              # lag 0
  "Isotherm_21_lag1",     # lag 1
  "Salinity_bottles",     # lag 0
  "u10_lag2",             # lag 2
  "tp_lag1",              # lag 1
  "AMO_lag3",             # lag 3
  "MEIv2_lag5"            # lag 3
)

# Store in list for looping
season_env_list <- list(
  "Upwelling" = GF_env_Upwelling,
  "Secondary Upwelling" = GF_env_SecondaryUpwelling,
  "Rainy" = GF_env_Rainy
)

# Run final models and store results
gf_results <- list()
importance_plots <- list()

for (s in names(season_env_list)) {
  cat("\n========================================\n")
  cat("Running final model for:", s, "\n")
  cat("========================================\n")
  
  env_factors <- season_env_list[[s]]
  
  # Extract data
  data <- extractMatrixFix_seasonal(env_factors, season_filter = s)
  envGF <- data[[1]]
  specGF <- data[[2]]
  
  # Filter genera
  mo5_specGF <- specGF %>% select_if(colSums(. > 0) > 10)
  colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))
  
  nSites <- dim(mo5_specGF)[1]
  nSpecs <- dim(mo5_specGF)[2]
  lev <- floor(log2(nSites * 0.368/2))
  
  cat("nSites:", nSites, "\n")
  cat("nSpecs:", nSpecs, "\n")
  cat("maxLevel:", lev, "\n")
  
  # Log transform
  log_specGF <- log1p(mo5_specGF * 100)
  
  # Run GF
  gf <- gradientForest(cbind(envGF, log_specGF),
                       predictor.vars = colnames(envGF),
                       response.vars = colnames(mo5_specGF),
                       ntree = 1500, transform = NULL, compact = TRUE,
                       trace = TRUE, nbin = 201, maxLevel = lev,
                       corr.threshold = 0.5)
  
  gf_results[[s]] <- gf
  
  # Print importance
  cat("\nImportance ranking:\n")
  print(importance(gf))
  
  # Create importance plot (matching your original style)
  imp.w <- importance(gf, "Weighted")
  o.w <- order(imp.w)
  
  # Clean up variable names for plotting
  var_names <- names(imp.w)
  var_names <- gsub("_merged", "", var_names)
  var_names <- gsub("_bottles", "", var_names)
  var_names <- gsub("_10m", "", var_names)
  var_names <- gsub("_lag", " lag", var_names)
  
  data_df <- data.frame(
    var = factor(var_names, levels = var_names[o.w]),
    importance = imp.w,
    row.names = NULL
  )
  
  title_expr <- expression(paste(R^2, " weighted importance"))
  
  p <- ggplot(data_df, aes(x = var, y = importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(data_df$var))) +
    xlab("Predictor Variable") +
    ylab(title_expr) +
    ggtitle(s) +
    theme_cowplot()
  
  importance_plots[[s]] <- p
  #print(p)
  
  # Save model
  season_label <- gsub(" ", "_", s)
  #saveRDS(gf, paste0("GF_final_model_", season_label, ".RDS"))
}

# Combined plot
combined_plot <- plot_grid(
  importance_plots[["Upwelling"]], 
  importance_plots[["Secondary Upwelling"]], 
  importance_plots[["Rainy"]],
  ncol = 3, labels = c("A", "B", "C")
)

print(combined_plot)
# ggsave("GF_seasonal_importance_comparison.pdf", combined_plot, width = 15, height = 6)
