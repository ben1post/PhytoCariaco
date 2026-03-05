library(tidyverse)
library(data.table)
library(vegan)
library(viridis)
library(cowplot)
library(gradientForest)
library(RColorBrewer)

# Read Phytoplankton Data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts <- phyto_counts[complete.cases(phyto_counts$counts),]
phyto_counts[phyto_counts$Genus=="Pseudo-nitzschia",]$Genus <- "Pseudo.nitzschia"

# Get environmental data
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

# Define variable sets (excluding evaporation "e")
GF_env_nolag <- c("u10", "NO3_merged", "PO4_merged", "SiO4_merged", 
                  "tp", "Salinity_bottles", "sst_10m", 
                  "Isotherm_21", "AMO", "MEIv2")

# Upwelling: optimal lags
GF_env_Upwelling <- c(
  "NO3_merged", "PO4_merged", "SiO4_merged",
  "sst_10m", "Isotherm_21", "Salinity_bottles",
  "u10", "tp_lag1", "AMO_lag3", "MEIv2_lag4"
)

# Secondary Upwelling: optimal lags
GF_env_SecondaryUpwelling <- c(
  "NO3_merged", "PO4_merged", "SiO4_merged_lag1",
  "sst_10m", "Isotherm_21", "Salinity_bottles_lag2",
  "u10", "tp_lag2", "AMO_lag2", "MEIv2_lag5"
)

# Rainy: optimal lags
GF_env_Rainy <- c(
  "NO3_merged", "PO4_merged", "SiO4_merged_lag2",
  "sst_10m", "Isotherm_21", "Salinity_bottles",
  "u10", "tp", "AMO_lag3", "MEIv2_lag5"
)

# Store configurations
season_config <- list(
  "Upwelling" = list(nolag = GF_env_nolag, lags = GF_env_Upwelling),
  "Secondary Upwelling" = list(nolag = GF_env_nolag, lags = GF_env_SecondaryUpwelling),
  "Rainy" = list(nolag = GF_env_nolag, lags = GF_env_Rainy)
)

# Function to run GF
run_seasonal_GF <- function(env_factors, season_filter) {
  data <- extractMatrixFix_seasonal(env_factors, season_filter)
  envGF <- data[[1]]
  specGF <- data[[2]]
  
  mo5_specGF <- specGF %>% select_if(colSums(. > 0) > 10)
  colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))
  
  nSites <- dim(mo5_specGF)[1]
  lev <- floor(log2(nSites * 0.368/2))
  
  log_specGF <- log1p(mo5_specGF) # * 100)
  
  gf <- gradientForest(cbind(envGF, log_specGF),
                       predictor.vars = colnames(envGF),
                       response.vars = colnames(mo5_specGF),
                       ntree = 1500, transform = NULL, compact = TRUE,
                       trace = TRUE, nbin = 201, maxLevel = lev,
                       corr.threshold = 0.5)
  
  return(list(gf = gf, specGF = mo5_specGF))
}

# --- Updated Variable Importance Function ---
create_var_imp_plot <- function(gf, title_text, x_limit = NULL, show_axis_labels = TRUE) {
  imp.w <- importance(gf, "Weighted")
  o.w <- order(imp.w, decreasing = TRUE)
  
  var_names <- names(imp.w)[o.w]
  var_names <- gsub("_merged", "", var_names)
  var_names <- gsub("_bottles", "", var_names)
  var_names <- gsub("_10m", "", var_names)
  var_names <- gsub("_lag", " L", var_names)
  
  data_df <- data.frame(
    var = factor(var_names, levels = rev(var_names)),
    importance = as.numeric(imp.w[o.w])
  )
  
  p <- ggplot(data_df, aes(x = var, y = importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("") +
    ylab(expression(paste(R^2, " weighted importance"))) +
    ggtitle(title_text) +
    theme_cowplot(font_size = 12) +
    theme(plot.title = element_text(size = 10))
  
  # Apply limits
  if (!is.null(x_limit)) {
    p <- p + ylim(0, x_limit) 
  }
  
  # Conditionally hide horizontal axis labels/title
  if (!show_axis_labels) {
    p <- p + theme(
      axis.title.x = element_blank()
    )
  }
  
  return(p)
}

# --- Updated Genus Importance Function ---
create_spec_imp_plot <- function(gf, gen_lookup, fg_col_scale, show_axis_labels = TRUE) {
  perf <- importance(gf, type = "Species")
  o.s <- order(perf)
  specnames <- names(perf[o.s])
  
  speccdata_df <- data.frame(
    Genus = specnames,
    importance = as.numeric(perf[o.s])
  )
  
  speccdata_df <- speccdata_df %>%
    left_join(gen_lookup %>% select(Genus, FuncGroup), by = "Genus") %>%
    mutate(Genus = factor(Genus, levels = specnames))
  
  p <- ggplot(speccdata_df, aes(x = Genus, y = importance, col = FuncGroup)) +
    geom_point(size = 2) +
    coord_flip() +
    ylab(expression(paste(R^2, " weighted importance"))) +
    xlab("") +
    ggtitle("Genus Importance (opt. lags)") +
    scale_colour_manual(values = fg_col_scale, name = "Functional Group") +
    theme_cowplot(font_size = 10) +
    theme(
      plot.title = element_text(size = 10),
      axis.text.y = element_text(size = 6),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "grey90"),
      legend.margin = margin(4, 4, 4, 4),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
  
  # Conditionally hide horizontal axis labels/title
  if (!show_axis_labels) {
    p <- p + theme(
      axis.title.x = element_blank()
    )
  }
  
  return(p)
}

# --- Color Scale Setup (Run this once before the loop) ---

# 1. Get unique functional groups per Genus
GenusFuncGroup_match <- phyto_counts %>% 
  group_by(Genus) %>% 
  summarize(FuncGroup = first(FuncGroup))

# 2. Define colors for Functional Groups
gens = unique(GenusFuncGroup_match$FuncGroup)
genscolscale = brewer.pal(length(gens), "Set1")
names(genscolscale) <- gens

# 3. Map these colors to the Genus level
GenusFuncGroup_match$color <- genscolscale[GenusFuncGroup_match$FuncGroup]

# 4. Create the final named vector for ggplot
speccolscale = as.character(GenusFuncGroup_match$color)
names(speccolscale) <- as.character(GenusFuncGroup_match$Genus)



season_plots_list <- list()
seasons <- c("Upwelling", "Secondary Upwelling", "Rainy")

# Limits from previous step
season_limits <- list(
  "Upwelling" = 0.04,
  "Secondary Upwelling" = 0.08,
  "Rainy" = 0.04
)

# Initialize list to store text summaries
summary_text_list <- list()

for (s in seasons) {
  cat("Processing season:", s, "\n")
  
  # --- Setup ---
  result_nolag <- run_seasonal_GF(season_config[[s]]$nolag, s)
  gf_nolag <- result_nolag$gf
  result_lags <- run_seasonal_GF(season_config[[s]]$lags, s)
  gf_lags <- result_lags$gf
  
  current_limit <- season_limits[[s]]
  
  # --- Determine if we show axis labels ---
  # Only show labels for the bottom-most plot (Rainy)
  show_labels <- (s == "Rainy")
  
  # --- Create plots with the new argument ---
  p_var_nolag <- create_var_imp_plot(gf_nolag, "Variable Importance (no lag)", 
                                     x_limit = current_limit, 
                                     show_axis_labels = show_labels)
  
  p_var_lags <- create_var_imp_plot(gf_lags, "Variable Importance (opt. lags)", 
                                    x_limit = current_limit, 
                                    show_axis_labels = show_labels)
  
  p_spec <- create_spec_imp_plot(gf_lags, GenusFuncGroup_match, genscolscale, 
                                 show_axis_labels = show_labels)
  
  # --- Combine Row ---
  combined_row <- plot_grid(p_var_nolag, p_var_lags, p_spec, 
                            ncol = 3, rel_widths = c(1, 1, 1.2), align = "h")
  
  # --- Add Title ---
  title <- ggdraw() + 
    draw_label(s, fontface = "bold", size = 16, x = 0.01, hjust = 0)
  
  season_panel <- plot_grid(title, combined_row, ncol = 1, rel_heights = c(0.1, 1))
  season_plots_list[[s]] <- season_panel
  
  
  # --- EXTRACT RESULTS SUMMARY (Final Version) ---
  
  # 1. Get Data Dimensions
  n_samples <- nrow(result_lags$specGF)      # Time points
  n_genus_input <- ncol(result_lags$specGF)  # Total species passed to GF
  n_preds   <- length(season_config[[s]]$lags)
  
  # 2. Get Species Performance (R2)
  spec_r2 <- importance(gf_lags, type = "Species")
  # Count species with positive variance explained
  n_genus_predicted <- sum(spec_r2 > 0) 
  
  # 3. Get Weighted Predictor Importance
  imp_vals <- importance(gf_lags, "Weighted")
  top5_indices <- order(imp_vals, decreasing = TRUE)[1:5]
  top5_names <- names(imp_vals)[top5_indices]
  top5_scores <- imp_vals[top5_indices]
  
  # 4. Format the text for predictors
  pred_summaries <- sapply(seq_along(top5_names), function(i) {
    name <- top5_names[i]
    score <- round(top5_scores[i], 5) 
    
    if (grepl("_lag", name)) {
      lag_num <- sub(".*_lag", "", name)
      base_var <- sub("_lag.*", "", name)
      base_var <- gsub("_merged", "", base_var)
      base_var <- gsub("_bottles", "", base_var)
      base_var <- gsub("_10m", "", base_var)
      return(paste0("   ", i, ". ", base_var, " (Lag: ", lag_num, " mo) | R2: ", score))
    } else {
      base_var <- name
      base_var <- gsub("_merged", "", base_var)
      base_var <- gsub("_bottles", "", base_var)
      base_var <- gsub("_10m", "", base_var)
      return(paste0("   ", i, ". ", base_var, " (No Lag)    | R2: ", score))
    }
  })
  
  # 5. Compile the text block
  summary_text_list[[s]] <- paste0(
    "SEASON: ", s, "\n",
    "--------------------------------------------------\n",
    "Dataset Stats:\n",
    "   - Samples (n):      ", n_samples, "\n",
    "   - Predictors:       ", n_preds, "\n",
    "   - Genus Input:      ", n_genus_input, "\n",
    "   - Genus Predicted:  ", n_genus_predicted, " (R2 > 0)\n",
    "--------------------------------------------------\n",
    "Top 5 Predictors (Weighted Importance):\n",
    paste(pred_summaries, collapse = "\n"),
    "\n\n"
  )
  
}

# --- Combine and Save ---
final_combined_figure <- plot_grid(
  plotlist = season_plots_list, 
  ncol = 1, 
  align = "v"
)
print(final_combined_figure)
# --- Save Final Plot ---
# Note: Height is increased to 24 (8 * 3) to accommodate three rows
#ggsave("plots/exports/GF_seasonal_ALL_COMBINED.pdf", 
#       final_combined_figure, width = 8.49, height = 10.57)



# --- PRINT TEXT SUMMARY FOR MANUSCRIPT ---
cat("\n\n##############################################\n")
cat("       SUMMARY OF RESULTS (Optimized Lags)    \n")
cat("##############################################\n\n")

for (season_txt in summary_text_list) {
  cat(season_txt)
}

