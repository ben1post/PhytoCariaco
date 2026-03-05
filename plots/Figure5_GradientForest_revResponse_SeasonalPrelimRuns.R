##### Seasonal Lag Analysis - Prelim Runs #####

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



# Run GF model - FIXED to match non-seasonal methodology
runGFmodel_seasonal <- function(env_factors, season_filter = NULL) {
  data <- extractMatrixFix_seasonal(env_factors, season_filter)
  
  envGF <- data[[1]]
  specGF <- data[[2]]
  
  mo5_specGF <- specGF %>% select_if(colSums(. > 0) > 5)
  colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))
  
  nSites <- dim(mo5_specGF)[1]
  lev <- floor(log2(nSites * 0.368/2))
  
  log_specGF <- log1p(mo5_specGF)
  
  gf <- gradientForest(cbind(envGF, log_specGF),
                       predictor.vars = colnames(envGF),
                       response.vars = colnames(mo5_specGF),
                       ntree = 1500, transform = NULL, compact = TRUE,
                       trace = TRUE, nbin = 201, maxLevel = lev,
                       corr.threshold = 0.5)
  return(gf)
}

# ------------------------------------------------------------------------------
# Function to run lag analysis for all seasons (returns variable_selection per season)
# ------------------------------------------------------------------------------
run_lag_analysis_all_seasons <- function(seasons = c("Upwelling", "Secondary Upwelling", "Rainy")) {
  
  results <- list()
  
  for (s in seasons) {
    cat("\nRunning lag analysis for:", s, "\n")
    
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
    Isotherm_21_lags <- safe_run(GF_env_factors_Isotherm, "Isotherm")
    Salinity_lags <- safe_run(GF_env_factors_Salinity, "Salinity")
    u10_lags <- safe_run(GF_env_factors_u10, "u10")
    tp_lags <- safe_run(GF_env_factors_tp, "tp")
    AMO_lags <- safe_run(GF_env_factors_AMO, "AMO")
    MEIv2_lags <- safe_run(GF_env_factors_MEIv2, "MEIv2")
    
    sq <- seq(17)
    
    variable_selection <- data.frame(
      lag = c(0:12, 24, 36, 48, 60),
      NO3 = NO3_lags[sq],
      PO4 = PO4_lags[sq],
      SiO4 = SiO4_lags[sq],
      Temperature = SST_lags[sq],
      Isotherm_21 = Isotherm_21_lags[sq],
      Salinity_bottles = Salinity_lags[sq],
      Wind_speed = u10_lags[sq],
      Precipitation = tp_lags[sq],
      AMO = AMO_lags[sq],
      MEIv2 = MEIv2_lags[sq]
    )
    
    results[[s]] <- variable_selection
  }
  
  return(results)
}

# ------------------------------------------------------------------------------
# 5. DEFINE LAG RANGES FOR EACH VARIABLE
# ------------------------------------------------------------------------------

# Local variables: short lags (0-3 months)
GF_env_factors_NO3 <- c("NO3_merged", "NO3_merged_lag1", "NO3_merged_lag2", "NO3_merged_lag3")
GF_env_factors_PO4 <- c("PO4_merged", "PO4_merged_lag1", "PO4_merged_lag2", "PO4_merged_lag3")
GF_env_factors_SiO4 <- c("SiO4_merged", "SiO4_merged_lag1", "SiO4_merged_lag2", "SiO4_merged_lag3")
GF_env_factors_SST <- c("sst_10m", "sst_10m_lag1", "sst_10m_lag2", "sst_10m_lag3")
GF_env_factors_Isotherm <- c("Isotherm_21", "Isotherm_21_lag1", "Isotherm_21_lag2", "Isotherm_21_lag3")
GF_env_factors_Salinity <- c("Salinity_bottles", "Salinity_bottles_lag1", "Salinity_bottles_lag2", "Salinity_bottles_lag3")

# Atmospheric variables: medium lags (0-6 months)
GF_env_factors_u10 <- c("u10", "u10_lag1", "u10_lag2", "u10_lag3", "u10_lag4", "u10_lag5", "u10_lag6")
GF_env_factors_tp <- c("tp", "tp_lag1", "tp_lag2", "tp_lag3", "tp_lag4", "tp_lag5", "tp_lag6")

# Climate indices: extended lags (0-12, 24, 36, 48, 60 months) - matching non-seasonal
GF_env_factors_AMO <- c("AMO", paste0("AMO_lag", c(1:12, 24, 36, 48, 60)))
GF_env_factors_MEIv2 <- c("MEIv2", paste0("MEIv2_lag", c(1:12, 24, 36, 48, 60)))


# ------------------------------------------------------------------------------
# GENERATE LATEX TABLES FOR SEASONAL ANALYSIS
# ------------------------------------------------------------------------------

library(kableExtra)
library(knitr)

options(scipen = 3)

# Function to bold the max value in a column
bold_max <- function(x) {
  cell_spec(round(x, 4), format = "latex", 
            bold = (x == max(x, na.rm = TRUE) & !is.na(x)))
}

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
  
  # Run models
  NO3_lags <- safe_run(GF_env_factors_NO3, "NO3")
  PO4_lags <- safe_run(GF_env_factors_PO4, "PO4")
  SiO4_lags <- safe_run(GF_env_factors_SiO4, "SiO4")
  SST_lags <- safe_run(GF_env_factors_SST, "SST")
  Isotherm_21_lags <- safe_run(GF_env_factors_Isotherm, "Isotherm")
  Salinity_lags <- safe_run(GF_env_factors_Salinity, "Salinity")
  u10_lags <- safe_run(GF_env_factors_u10, "u10")
  tp_lags <- safe_run(GF_env_factors_tp, "tp")
  AMO_lags <- safe_run(GF_env_factors_AMO, "AMO")
  MEIv2_lags <- safe_run(GF_env_factors_MEIv2, "MEIv2")
  
  # --------------------------------------------------------------------------
  # BUILD COMPARISON TABLE
  # --------------------------------------------------------------------------
  
  sq <- seq(17)
  
  variable_selection <- data.frame(
    lag = c(0:12, 24, 36, 48, 60),
    NO3 = NO3_lags[sq],
    PO4 = PO4_lags[sq],
    SiO4 = SiO4_lags[sq],
    Temperature = SST_lags[sq],
    Isotherm_21 = Isotherm_21_lags[sq],
    Salinity_bottles = Salinity_lags[sq],
    Wind_speed = u10_lags[sq],
    Precipitation = tp_lags[sq],
    AMO = AMO_lags[sq],
    MEIv2 = MEIv2_lags[sq]
  )
  
  print(variable_selection)
  
  # --------------------------------------------------------------------------
  # GENERATE LATEX TABLE
  # --------------------------------------------------------------------------
  
  season_label <- tolower(gsub(" ", "_", s))
  
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
  
  col_names <- c("Lag", "NO\\textsubscript{3}", "PO\\textsubscript{4}", 
                 "SiO\\textsubscript{4}", "SST", 
                 "\\qty{21}{\\celsius} Isotherm", "Salinity",
                 "Wind speed", "Precipitation", "AMO", "MEI v.2")
  
  cat("\n\n% === LATEX TABLE:", s, "===\n")
  variable_selection_formatted %>%
    kable(format = "latex", booktabs = TRUE, escape = FALSE,
          col.names = col_names,
          caption = paste0("$R^2$ weighted importance for predictor variables during the ", 
                           s, " season. In-situ variables tested with lags 0--3 months, ",
                           "atmospheric variables 0--6 months, and climate indices 0--12, 24, 36, 48, 60 months. ",
                           "Maximum importance per variable highlighted in bold."),
          label = paste0("tab:sup:GFlagtests_", season_label)) %>%
    kable_styling(latex_options = "scale_down") %>%
    print()
}








# ##############################################################################
# AVERAGED LAG ANALYSIS (5 RUNS) AND LATEX TABLE GENERATION
# ##############################################################################

n_runs <- 5
seasons <- c("Upwelling", "Secondary Upwelling", "Rainy")

# Storage for all runs
all_results <- vector("list", n_runs)

for (run in 1:n_runs) {
  cat("\n\n############################################################\n")
  cat("RUN", run, "of", n_runs, "\n")
  cat("############################################################\n")
  
  all_results[[run]] <- run_lag_analysis_all_seasons(seasons)
}

# ------------------------------------------------------------------------------
# Average results across runs
# ------------------------------------------------------------------------------
averaged_results <- list()

for (s in seasons) {
  run_matrices <- lapply(all_results, function(x) {
    df <- x[[s]]
    as.matrix(df[, -1])  # exclude lag column
  })
  
  avg_matrix <- Reduce(`+`, run_matrices) / n_runs
  
  averaged_results[[s]] <- data.frame(
    lag = c(0:12, 24, 36, 48, 60),
    avg_matrix
  )
  colnames(averaged_results[[s]]) <- colnames(all_results[[1]][[s]])
}
# ------------------------------------------------------------------------------
# Generate formatted dataframe for LaTeX table
# ------------------------------------------------------------------------------
format_variable_cell <- function(importance_vector, lag_vector, threshold = 0.15) {
  valid_idx <- !is.na(importance_vector)
  if (sum(valid_idx) == 0) return("")
  
  imp <- importance_vector[valid_idx]
  lags <- lag_vector[valid_idx]
  
  # Top 3 by importance (descending)
  top_order <- order(imp, decreasing = TRUE)
  top_n <- min(3, length(top_order))
  top_idx <- top_order[1:top_n]
  
  top_lags <- lags[top_idx]
  top_imps <- imp[top_idx]
  
  # Selected lag: shortest within threshold of max
  max_imp <- max(imp)
  threshold_value <- max_imp * (1 - threshold)
  within_threshold <- which(imp >= threshold_value)
  selected_lag <- min(lags[within_threshold])
  
  entries <- sapply(1:top_n, function(i) {
    lag_val <- top_lags[i]
    imp_val <- sprintf("%.3f", top_imps[i])
    
    if (lag_val == selected_lag) {
      paste0("\\textbf{", lag_val, "} (", imp_val, ")")
    } else {
      paste0(lag_val, " (", imp_val, ")")
    }
  })
  
  paste(entries, collapse = "; ")
}

variable_meta <- list(
  list(col = "NO3", display = "NO\\textsubscript{3}", category = "In-situ", lags = 0:3),
  list(col = "PO4", display = "PO\\textsubscript{4}", category = "In-situ", lags = 0:3),
  list(col = "SiO4", display = "SiO\\textsubscript{4}", category = "In-situ", lags = 0:3),
  list(col = "Temperature", display = "SST", category = "In-situ", lags = 0:3),
  list(col = "Isotherm_21", display = "\\qty{21}{\\celsius} Isotherm", category = "In-situ", lags = 0:3),
  list(col = "Salinity_bottles", display = "Salinity", category = "In-situ", lags = 0:3),
  list(col = "Wind_speed", display = "Wind speed", category = "Atmospheric", lags = 0:6),
  list(col = "Precipitation", display = "Precipitation", category = "Atmospheric", lags = 0:6),
  list(col = "AMO", display = "AMO", category = "Climate", lags = c(0:12, 24, 36, 48, 60)),
  list(col = "MEIv2", display = "MEI v.2", category = "Climate", lags = c(0:12, 24, 36, 48, 60))
)

seasons <- c("Upwelling", "Secondary Upwelling", "Rainy")

# Build the summary dataframe
lag_summary_table <- do.call(rbind, lapply(variable_meta, function(v) {
  cells <- sapply(seasons, function(s) {
    df <- averaged_results[[s]]
    lag_idx <- which(df$lag %in% v$lags)
    format_variable_cell(df[[v$col]][lag_idx], df$lag[lag_idx], threshold = 0.20)
  })
  
  data.frame(
    Variable = v$display,
    Category = v$category,
    Upwelling = cells[1],
    Secondary_Upwelling = cells[2],
    Rainy = cells[3],
    stringsAsFactors = FALSE
  )
}))


lag_summary_table <- data.frame(
  Variable = c(
    "NO\\textsubscript{3}",
    "PO\\textsubscript{4}",
    "SiO\\textsubscript{4}",
    "SST",
    "\\qty{21}{\\celsius} Isotherm",
    "Salinity",
    "Wind speed",
    "Precipitation",
    "AMO",
    "MEI v.2"
  ),
  Category = c(
    rep("In-situ", 6),
    rep("Atmospheric", 2),
    rep("Climate", 2)
  ),
  Upwelling = c(
    "\\textbf{0} (0.060); 1 (0.036); 2 (0.020)",
    "\\textbf{0} (0.023); 1 (0.008); 2 (0.007)",
    "\\textbf{0} (0.026); 1 (0.008); 3 (0.006)",
    "\\textbf{0} (0.101); 3 (0.014); 1 (0.007)",
    "\\textbf{0} (0.045); 1 (0.011); 2 (0.007)",
    "\\textbf{0} (0.031); 3 (0.010); 2 (0.008)",
    "5 (0.007); 6 (0.007); 4 (0.006)",
    "\\textbf{1} (0.019); 0 (0.012); 2 (0.009)",
    "60 (0.028); \\textbf{24} (0.023); 36 (0.021)",
    "\\textbf{48} (0.036); 36 (0.018); 60 (0.018)"
  ),
  Secondary_Upwelling = c(
    "\\textbf{0} (0.034); 3 (0.031); 2 (0.029)",
    "3 (0.023); \\textbf{0} (0.021); 2 (0.020)",
    "\\textbf{0} (0.032); 1 (0.028); 2 (0.022)",
    "\\textbf{0} (0.032); 2 (0.025); 1 (0.014)",
    "\\textbf{1} (0.024); 3 (0.014); 0 (0.010)",
    "\\textbf{2} (0.073); 3 (0.049); 1 (0.020)",
    "\\textbf{6} (0.027); 2 (0.010); 5 (0.006)",
    "\\textbf{0} (0.018); 2 (0.012); 1 (0.008)",
    "\\textbf{24} (0.034); 60 (0.032); 36 (0.020)",
    "\\textbf{60} (0.037); 24 (0.018); 0 (0.015)"
  ),
  Rainy = c(
    "3 (0.017); \\textbf{0} (0.015); 2 (0.006)",
    "\\textbf{0} (0.039); 3 (0.027); 1 (0.014)",
    "\\textbf{3} (0.043); 1 (0.033); 2 (0.032)",
    "\\textbf{0} (0.052); 1 (0.032); 3 (0.010)",
    "\\textbf{0} (0.032); 2 (0.010); 3 (0.007)",
    "\\textbf{0} (0.022); 3 (0.010); 2 (0.010)",
    "\\textbf{0} (0.021); 5 (0.019); 2 (0.009)",
    "\\textbf{3} (0.017); 0 (0.012); 6 (0.006)",
    "48 (0.028); \\textbf{36} (0.025); 12 (0.013)",
    "48 (0.027); \\textbf{36} (0.024); 60 (0.020)"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Generate LaTeX table using kable
# ------------------------------------------------------------------------------
library(kableExtra)
library(knitr)

lag_summary_table %>%
  select(-Category) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE,
        row.names = FALSE,
        col.names = c("Variable", "Upwelling", "Secondary Upwelling", "Rainy"),
        caption = "Mean $R^2$ weighted importance for predictor variables across three seasons, used to select optimal time lags. Values shown as: lag in months (mean importance). The analysis was run five times per season due to relatively small sample sizes, and mean importance values are reported. For in-situ variables, lags of 0--3 months were tested; for atmospheric variables, 0--6 months; for climate indices, 0--12, 24, 36, 48, and 60 months. The shortest lag within 20\\% of the maximum importance is highlighted in bold and selected for final models.",
        label = "sup:GFlagtests_seasonal") %>%
  kable_styling(latex_options = "scale_down") %>%
  pack_rows("In-situ variables (lags 0--3 months)", 1, 6, italic = TRUE, bold = FALSE) %>%
  pack_rows("Atmospheric variables (lags 0--6 months)", 7, 8, italic = TRUE, bold = FALSE) %>%
  pack_rows("Climate indices (lags 0--12, 24, 36, 48, 60 months)", 9, 10, italic = TRUE, bold = FALSE)











