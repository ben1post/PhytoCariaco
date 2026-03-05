library(tidyverse)
library(cowplot)
library(scales)


library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)
library(lubridate)

# --- 1. DATA PREP (Original) ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds") 

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep=""), format="%m-%Y-%d"),
    year = as.numeric(format(date, "%Y")),
    month_num = as.numeric(format(date, "%m")),
    month_label = factor(month.abb[month_num], levels = month.abb)
  )

# --- 2. DATA PREP (Cluster Groups) ---
# Update group names to match your requested legend labels exactly
GroupCluster1 <- CARIACO %>%
  filter((date >= as.Date("1996-01-01") & date <= as.Date("2003-12-31")) | 
           (date >= as.Date("2014-01-01") & date <= as.Date("2016-12-31"))) %>%
  mutate(group = "Cluster 1")

GroupCluster2 <- CARIACO %>%
  filter((date >= as.Date("2004-01-01") & date <= as.Date("2013-12-31")) | 
           (date >= as.Date("2017-06-01") & date <= as.Date("2017-12-31"))) %>%
  mutate(group = "Cluster 2")

ENV_DATA_groups <- rbind(GroupCluster1, GroupCluster2)

# Define the exact color palette requested
cluster_cols <- c("Cluster 1" = "blue", 
                  "Cluster 2" = "red")

library(tidyverse)
library(cowplot)
library(scales)

# --- 1. DATA PREP ---

# Helper variables to summarize
target_vars <- c("Isotherm_21", "NO3_merged", "Chlorophyll")

# A. Calculate Total Mean & SD (Black Line / Grey Ribbon)
stats_total <- CARIACO %>%
  group_by(month_label) %>%
  summarise(across(all_of(target_vars), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  mutate(group = "Total Mean")

# B. Calculate Cluster Means & SD (Blue/Red Lines & Ribbons)
stats_cluster <- ENV_DATA_groups %>%
  group_by(month_label, group) %>%
  summarise(across(all_of(target_vars), 
                   list(mean = ~mean(., na.rm = TRUE), 
                        sd = ~sd(., na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

# --- 2. PLOTTING FUNCTION ---

create_season_plot <- function(var_root, y_lab, plot_title = NULL, show_legend = FALSE, show_x = FALSE, invert_y = FALSE) {
  
  # Construct column names dynamically (e.g., "Chlorophyll_mean", "Chlorophyll_sd")
  mean_col <- sym(paste0(var_root, "_mean"))
  sd_col   <- sym(paste0(var_root, "_sd"))
  
  p <- ggplot() +
    
    # 1. Ribbon for Total Mean (Grey background context)
    geom_ribbon(data = stats_total, 
                aes(x = month_label, ymin = !!mean_col - !!sd_col, ymax = !!mean_col + !!sd_col, group = 1),
                fill = "grey40", alpha = 0.15) +
    
    # 2. Ribbons for Clusters (Blue/Red Variance)
    geom_ribbon(data = stats_cluster, 
                aes(x = month_label, ymin = !!mean_col - !!sd_col, ymax = !!mean_col + !!sd_col, 
                    fill = group, group = group),
                alpha = 0.15) +
    
    # 3. Lines for Clusters
    geom_line(data = stats_cluster, 
              aes(x = month_label, y = !!mean_col, color = group, group = group), 
              size = 1) +
    geom_point(data = stats_cluster, 
               aes(x = month_label, y = !!mean_col, color = group), 
               size = 1.5) +
    
    # 4. Line for Total Mean (Black)
    geom_line(data = stats_total, 
              aes(x = month_label, y = !!mean_col, group = 1), 
              color = "black", size = 1.2, linetype = "solid") + # Changed to solid for clarity over ribbons
    
    # Styling
    scale_color_manual(values = cluster_cols) +
    scale_fill_manual(values = cluster_cols) +
    theme_cowplot() +
    labs(y = y_lab, title = plot_title) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(b = 5, t = 5)
    )
  
  # Conditional: Invert Y Axis (Good for Depth)
  if (invert_y) {
    p <- p + scale_y_reverse()
  }
  
  # Conditional: Show X Axis Labels
  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  # Conditional: Show Legend
  if (show_legend) {
    p <- p + theme(legend.position = "top", legend.justification = "center", legend.title = element_blank())
  }
  
  return(p)
}

# --- 3. GENERATE PANELS ---

# Panel 1: Isotherm 21 (Inverted Y-axis)
p_iso <- create_season_plot("Isotherm_21", "Iso 21 Depth (m)", 
                            plot_title = "Seasonal Climatology (+/- SD)",
                            show_legend = TRUE, 
                            invert_y = TRUE) # <--- Invert Triggered Here

# Panel 2: Nitrate
p_no3 <- create_season_plot("NO3_merged", expression(NO[3]~(mu*M)), 
                            show_legend = FALSE)

# Panel 3: Chlorophyll
p_chl <- create_season_plot("Chlorophyll", expression(Chl~a~(mg~m^{-3})), 
                            show_legend = FALSE, show_x = TRUE) 

# --- 4. ALIGN AND DISPLAY ---
aligned_plots <- align_plots(p_iso, p_no3, p_chl, align = "v", axis = "l")

final_stack <- plot_grid(
  aligned_plots[[1]],
  aligned_plots[[2]],
  aligned_plots[[3]],
  ncol = 1,
  rel_heights = c(1.3, 1, 1)
)

print(final_stack)
