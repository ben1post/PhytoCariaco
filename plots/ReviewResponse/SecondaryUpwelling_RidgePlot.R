library(tidyverse)
library(cowplot)
library(ggpubr)
library(egg)

# --- 1. LOAD & PREP DATA ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

# Prepare variables
CARIACO$u10_negative <- -CARIACO$u10

# --- 2. DEFINE SEASONS ---
get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  factor(season_lookup[months], levels = c("Upwelling", "Secondary Upwelling", "Rainy"))
}

# --- 3. ASSIGN CLUSTERS & SEASONS ---
CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    season = get_season_3cat(date),
    cluster = case_when(
      (date >= as.Date("1996-01-01") & date <= as.Date("2003-12-31")) |
        (date >= as.Date("2014-01-01") & date <= as.Date("2016-12-31")) ~ "Cluster 1",
      (date >= as.Date("2004-01-01") & date <= as.Date("2013-12-31")) |
        (date >= as.Date("2017-06-01") & date <= as.Date("2017-12-31")) ~ "Cluster 2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cluster))

# --- 4. DEFINE COLORS ---
cluster_cols <- c("Cluster 1" = "red", "Cluster 2" = "blue")

# --- 5. UPDATED PLOTTING FUNCTION ---
make_density_plot <- function(data, x_var, x_lab, log_scale = FALSE, reverse_scale = FALSE, show_legend = FALSE) {
  
  # 1. Calculate Medians manually per Season AND Cluster
  # We use .data[[x_var]] to dynamically reference the column name passed as a string
  medians <- data %>%
    group_by(season, cluster) %>%
    summarise(med_val = median(.data[[x_var]], na.rm = TRUE), .groups = "drop")
  
  # 2. Base ggdensity (without built-in median)
  p <- ggdensity(
    data, 
    x = x_var,
    rug = TRUE,               # Keep the rug
    color = "cluster", 
    fill = "cluster",
    palette = cluster_cols,
    alpha = 0.5,
    bi = 1
  ) 
  
  # 3. Add the manually calculated medians
  # We use geom_vline with the 'medians' dataframe. 
  # Since 'medians' has a 'season' column, facet_grid will automatically map the right line to the right panel.
  p <- p + geom_vline(data = medians, aes(xintercept = med_val, color = cluster), 
                      linetype = "dashed", size = 0.4)
  
  # 4. Facet & Theme
  p <- p + 
    facet_grid(season ~ ., switch = "y") + 
    xlab(x_lab) +
    theme_cowplot() +
    theme(
      legend.position = if(show_legend) "right" else "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1)
    )
  
  # 5. Apply Scales
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  
  if (reverse_scale) {
    p <- p + scale_x_reverse()
  }
  
  return(p)
}

# --- 6. GENERATE PANELS ---

# 1. Wind Speed
p_wind <- make_density_plot(CARIACO, "u10_negative", "Wind Speed (m/s)")

# 2. Total Precipitation
p_tp <- make_density_plot(CARIACO, "tp", "Total Precip. (m)") +
  theme(strip.text.y.left = element_blank())

# 3. Iso 21 Depth (Reversed)
p_iso <- make_density_plot(CARIACO, "Isotherm_21", "Iso 21 Depth (m)", reverse_scale = TRUE) +
  theme(strip.text.y.left = element_blank())

# 4. Nitrate
p_no3 <- make_density_plot(CARIACO, "NO3_merged", expression(NO[3]~(mu*M))) +
  theme(strip.text.y.left = element_blank())

# 5. Chlorophyll (Log Scale + Legend)
p_pp <- make_density_plot(CARIACO, "PrimaryProductivity","PP", log_scale = TRUE) +
  theme(strip.text.y.left = element_blank())

# 5. Chlorophyll (Log Scale + Legend)
p_chl <- make_density_plot(CARIACO, "Chlorophyll", expression(Chl~a~(mg~m^{-3})), log_scale = TRUE, show_legend = TRUE) +
  theme(strip.text.y.left = element_blank())

# --- 7. ARRANGE ---
final_plot <- egg::ggarrange(p_wind, p_tp, p_iso, p_pp, p_chl, ncol = 5)
print(final_plot)







####### Cell Counts FG #######



# 1. Diatoms (Shows Season Labels)
p_diat <- make_density_plot(CARIACO, "Abundance_Diatom", expression(Diatoms~(Cells~ml^{-1})), log_scale = TRUE)

# 2. Haptophytes
p_hapt <- make_density_plot(CARIACO, "Abundance_Hapto", expression(Hapto~(Cells~ml^{-1})), log_scale = TRUE) +
  theme(strip.text.y.left = element_blank())

# 3. Dinoflagellates
p_dino <- make_density_plot(CARIACO, "Abundance_Dino", expression(Dino~(Cells~ml^{-1})), log_scale = TRUE) +
  theme(strip.text.y.left = element_blank())

# 4. Cyanobacteria
p_cyan <- make_density_plot(CARIACO, "Abundance_Cyano", expression(Cyano~(Cells~ml^{-1})), log_scale = TRUE) +
  theme(strip.text.y.left = element_blank())

# 5. Nanoflagellates (Shows Legend)
p_nano <- make_density_plot(CARIACO, "Abundance_Nanoflagellate", expression(Nano~(Cells~ml^{-1})), log_scale = TRUE, show_legend = TRUE) +
  theme(strip.text.y.left = element_blank())

# --- 7. ARRANGE ---
final_bio_plot <- egg::ggarrange(p_diat, p_hapt, p_dino, p_cyan, p_nano, ncol = 5)
print(final_bio_plot)
