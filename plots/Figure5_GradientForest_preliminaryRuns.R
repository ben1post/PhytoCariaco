# ==============================================================================
# Gradient Forest Lag Selection Analysis
# ==============================================================================
#
# Purpose: Determine optimal time lags for environmental predictors by running
#          separate GF models for each variable across a range of lags and
#          comparing weighted R² importance values.
#
# Output:  variable_selection dataframe comparing importance across lags
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(viridis)
library(cowplot)
library(gradientForest)

# ------------------------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------------------------

phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts), ]
phyto_counts[phyto_counts$Genus == "Pseudo-nitzschia", ]$Genus <- "Pseudo.nitzschia"

CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# ------------------------------------------------------------------------------
# 3. DATA EXTRACTION FUNCTION
# ------------------------------------------------------------------------------

#' Extract and align phytoplankton and environmental data for GF
extract_gf_matrices <- function(env_factors) {
  
  
  ds_genus <- phyto_counts %>%
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
    group_by(Genus, date) %>%
    summarise(Total = sum(counts), .groups = "drop") %>%
    arrange(desc(date))
  
  genus_wide <- pivot_wider(ds_genus,
                            names_from = Genus,
                            values_from = Total,
                            values_fn = sum,
                            values_fill = 0)
  genus_wide$time_month <- format(genus_wide$date, format = "%m-%Y")
  genus_wide <- genus_wide %>% select(-date)
  
  env_subset <- CARIACO %>% select(time_month, all_of(env_factors))
  
  combined <- list(env_subset, genus_wide) %>%
    reduce(full_join, by = "time_month") %>%
    select(-time_month)
  
  complete_data <- combined[complete.cases(combined), ]
  
  env_matrix <- complete_data %>% select(all_of(env_factors))
  genus_matrix <- complete_data %>% select(-all_of(env_factors))
  
  return(list(env_matrix, genus_matrix))
}

# ------------------------------------------------------------------------------
# 4. GRADIENT FOREST MODEL FUNCTION
# ------------------------------------------------------------------------------

#' Fit GF model for a set of lagged predictors
run_gf_model <- function(env_factors) {
  
  data <- extract_gf_matrices(env_factors)
  env_matrix <- data[[1]]
  genus_matrix <- data[[2]]
  
  # Filter rare genera
  genus_filtered <- genus_matrix %>% select_if(colSums(. > 0) > 10)
  colnames(genus_filtered) <- make.names(colnames(genus_filtered))
  
  # Calculate optimal tree depth
  n_sites <- nrow(genus_filtered)
  max_level <- floor(log2(n_sites * 0.368 / 2))
  
  # Log-transform counts
  genus_log <- log1p(genus_filtered)
  
  # Fit model
  gf <- gradientForest(cbind(env_matrix, genus_log),
                       predictor.vars = colnames(env_matrix),
                       response.vars = colnames(genus_filtered),
                       ntree = 1500,
                       transform = NULL,
                       compact = TRUE,
                       trace = TRUE,
                       nbin = 201,
                       maxLevel = max_level,
                       corr.threshold = 0.5)
  return(gf)
}

# ------------------------------------------------------------------------------
# 5. DEFINE LAG RANGES FOR EACH VARIABLE
# ------------------------------------------------------------------------------

# Local variables: short lags (0-3 months)
env_factors_NO3 <- c("NO3_merged", "NO3_merged_lag1", "NO3_merged_lag2", "NO3_merged_lag3")
env_factors_PO4 <- c("PO4_merged", "PO4_merged_lag1", "PO4_merged_lag2", "PO4_merged_lag3")
env_factors_SiO4 <- c("SiO4_merged", "SiO4_merged_lag1", "SiO4_merged_lag2", "SiO4_merged_lag3")
env_factors_SST <- c("sst_10m", "sst_10m_lag1", "sst_10m_lag2", "sst_10m_lag3")
env_factors_Isotherm <- c("Isotherm_21", "Isotherm_21_lag1", "Isotherm_21_lag2", "Isotherm_21_lag3")
env_factors_Salinity <- c("Salinity_bottles", "Salinity_bottles_lag1", "Salinity_bottles_lag2", "Salinity_bottles_lag3")

# Atmospheric variables: medium lags (0-6 months)
env_factors_u10 <- c("u10", "u10_lag1", "u10_lag2", "u10_lag3", "u10_lag4", "u10_lag5", "u10_lag6")
env_factors_tp <- c("tp", "tp_lag1", "tp_lag2", "tp_lag3", "tp_lag4", "tp_lag5", "tp_lag6")

# Climate indices: extended lags (0-60 months)
env_factors_AMO <- c("AMO", paste0("AMO_lag", c(1:12, 24, 36, 48, 60)))
env_factors_MEIv2 <- c("MEIv2", paste0("MEIv2_lag", c(1:12, 24, 36, 48, 60)))

# ------------------------------------------------------------------------------
# 6. RUN MODELS FOR EACH VARIABLE
# ------------------------------------------------------------------------------

cat("Running GF models for lag selection...\n\n")

gf_NO3 <- run_gf_model(env_factors_NO3)
gf_PO4 <- run_gf_model(env_factors_PO4)
gf_SiO4 <- run_gf_model(env_factors_SiO4)
gf_SST <- run_gf_model(env_factors_SST)
gf_Isotherm <- run_gf_model(env_factors_Isotherm)
gf_Salinity <- run_gf_model(env_factors_Salinity)
gf_u10 <- run_gf_model(env_factors_u10)
gf_tp <- run_gf_model(env_factors_tp)
gf_AMO <- run_gf_model(env_factors_AMO)
gf_MEIv2 <- run_gf_model(env_factors_MEIv2)

# ------------------------------------------------------------------------------
# 7. EXTRACT IMPORTANCE VALUES
# ------------------------------------------------------------------------------

# Extract importance in original order (not sorted by value)
NO3_lags <- as.numeric(importance(gf_NO3, sort = FALSE))
PO4_lags <- as.numeric(importance(gf_PO4, sort = FALSE))
SiO4_lags <- as.numeric(importance(gf_SiO4, sort = FALSE))
Temperature_lags <- as.numeric(importance(gf_SST, sort = FALSE))
Isotherm_21_lags <- as.numeric(importance(gf_Isotherm, sort = FALSE))
Salinity_lags <- as.numeric(importance(gf_Salinity, sort = FALSE))
u10_lags <- as.numeric(importance(gf_u10, sort = FALSE))
tp_lags <- as.numeric(importance(gf_tp, sort = FALSE))
AMO_lags <- as.numeric(importance(gf_AMO, sort = FALSE))
MEIv2_lags <- as.numeric(importance(gf_MEIv2, sort = FALSE))

# ------------------------------------------------------------------------------
# 8. BUILD COMPARISON TABLE
# ------------------------------------------------------------------------------

# Index for 17 lag positions (0-12, 24, 36, 48, 60)
sq <- seq(17)

variable_selection <- data.frame(
  lag = c(0:12, 24, 36, 48, 60),
  NO3 = NO3_lags[sq],
  PO4 = PO4_lags[sq],
  SiO4 = SiO4_lags[sq],
  Temperature = Temperature_lags[sq],
  Isotherm_21 = Isotherm_21_lags[sq],
  Salinity_bottles = Salinity_lags[sq],
  Wind_speed = u10_lags[sq],
  Precipitation = tp_lags[sq],
  AMO = AMO_lags[sq],
  MEIv2 = MEIv2_lags[sq]
)

print(variable_selection)


# ------------------------------------------------------------------------------
# 8. GENERATE LATEX TABLE FOR MANUSCRIPT
# ------------------------------------------------------------------------------


library(kableExtra)
library(knitr)

options(scipen=3)

# Function to bold the max value in a column
bold_max <- function(x) {
  cell_spec(round(x, 3), format = "latex", 
            bold = (x == max(x, na.rm = TRUE) & !is.na(x)))
}

# Reorder columns and rename
variable_selection_formatted <- variable_selection %>%
  select(lag, NO3, PO4, SiO4, 
         SST = Temperature, 
         Isotherm_21,
         Salinity = Salinity_bottles,
         Wind_speed,
         Precipitation,
         AMO, 
         MEIv2) %>%
  mutate(across(-lag, bold_max))

# Custom column names (with LaTeX commands) - matching the order above
col_names <- c("lag", "NO3", "PO4", "SiO4", "SST", 
               "\\qty{21}{\\celsius} Isotherm", "Salinity", 
               "Wind speed", "Precipitation", "AMO", "MEI v.2")

# Generate table
variable_selection_formatted %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE,
        col.names = col_names,
        caption = "Resulting $R^2$ weighted importance for running predictor variables together with time-lags to check if they provide a better fit to the data. Each column represents a single run, with resulting importance values for time series and monthly lags. The results were stable, so that repeated runs confirmed the chosen time lag. For in-situ variables, we tested all measurements up to a lag of 3 months. For climate variables we tested up to a lag of 6 months and added a lag of 12 months for the climate indices. Time lag with maximum importance per variable are highlighted in bold and chosen for the final model run.",
        label = "tab:sup:GFlagtests") %>%
  kable_styling(latex_options = "scale_down")



