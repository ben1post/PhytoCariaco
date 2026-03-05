library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)


# read combined dataset:
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")
CARIACO$year <- format(CARIACO$date, format="%Y")

# CARIACO$PrimaryProductivity
# units:milligrams Carbon/meter^3/hour (mgC/m^3/hr)

pp_yearly <- CARIACO %>%
  group_by(year) %>%
  filter(n() > 10) %>%
  summarize(mean_pp = mean(PrimaryProductivity, na.rm = TRUE))

pp_yearly$mean_pp / 1000 * 100 * 24*365




library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)
library(lubridate)

# --- 1. DATA PREP ---
# (Assuming 'CARIACO' is already loaded as per your snippet)

# Ensure date is parsed correctly and extract Month/Year factors
# We use factor() for Month to ensure it plots Jan->Dec, not alphabetical
CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep=""), format="%m-%Y-%d"),
    # Ensure year is numeric for continuous axis plotting
    year = as.numeric(format(date, "%Y")),
    month_num = as.numeric(format(date, "%m")),
    # Create ordered factor for month labels
    month_label = factor(month.abb[month_num], levels = month.abb)
  )

# --- 2. CALCULATE AVERAGE SEASONALITY ---
# This creates the "average year" to clearly show the two peaks
seasonal_clim <- CARIACO %>%
  group_by(month_label) %>%
  summarise(
    avg_chl = mean(PrimaryProductivity, na.rm = TRUE),
    # Calculate Standard Deviation for the ribbon
    sd_chl = sd(PrimaryProductivity, na.rm = TRUE)
  )

# --- 3. PLOT A: The Average Seasonal Cycle ---
# This proves your point about the Winter vs Summer peaks
p_season <- ggplot(seasonal_clim, aes(x = month_label, y = avg_chl)) +
  # Add a ribbon to show variability (+/- 1 SD)
  geom_ribbon(aes(ymin = pmax(0, avg_chl - sd_chl), # Ensure ribbon doesn't go below 0
                  ymax = avg_chl + sd_chl, group=1), 
              fill = "#2ca25f", alpha = 0.2) +
  geom_line(group = 1, color = "#2ca25f", size = 1.2) +
  geom_point(color = "#2ca25f", size = 3) +
  # Optional conceptual labels for the peaks
  #annotate("text", x = "Jan", y = max(seasonal_clim$avg_chl) * 1.05, 
  #         label = "Winter\nPeak", vjust = 0, size=3.5, fontface="bold", color="#006d2c") +
  #annotate("text", x = "Jul", y = seasonal_clim$avg_chl[seasonal_clim$month_label=="Jul"]*1.1, 
  #         label = "Summer\nPeak", vjust = 0, size=3.5, fontface="bold", color="#006d2c") +
  theme_cowplot() +
  labs(x = "Month", y = "Avg. Chlorophyll\n(mg m-3)", 
       title = "A. Average Seasonal Cycle (Climatology)") +
  theme(axis.text.x = element_blank(), # Hide X labels as they are in the plot below
        axis.title.x = element_blank(),
        plot.margin = margin(b = 2, unit = "pt")) # tighten margin between plots

# --- 4. PLOT B: The Time Series Heatmap (FLIPPED AXES) ---
# X = Month, Y = Year
p_heatmap <- ggplot(CARIACO, aes(x = month_label, y = year, fill = PrimaryProductivity)) +
  geom_tile(color = "white", size=0.2) + # Thin white borders for clarity
  # Use a nice green palette representing algae/chlorophyll. 
  # direction = 1 makes darker colors higher values.
  scale_fill_distiller(palette = "Greens", direction = 1, na.value = "grey95",
                       name = "Chl-a") +
  # Y Axis setup: Reverse trans puts newest years at the top
  scale_y_continuous(expand = c(0, 0), breaks = pretty_breaks(n = 10), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  theme_cowplot() +
  labs(x = "Month", y = "Year", title = "B. Full Time Series Overview") +
  theme(legend.position = "right",
        # Rotate x-axis labels for readability
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 2, unit = "pt"))

# --- 5. COMBINE AND ALIGN ---
# Aligning vertically ensures the month columns match perfectly
final_plot <- plot_grid(p_season, p_heatmap, ncol = 1, 
                        rel_heights = c(0.35, 0.65), align = "v", axis = "lr")

# Display
print(final_plot)




library(tidyverse)
library(cowplot)

# --- 1. PREPARE DATA ---
# (Assuming 'CARIACO' is already loaded)

# Calculate the correlation coefficient (Pearson) ignoring missing values
# We store it to print on the plot later
cor_res <- cor.test(CARIACO$Chlorophyll, CARIACO$PrimaryProductivity, 
                    method = "pearson", use = "complete.obs")
r_value <- round(cor_res$estimate, 2)
p_value <- format.pval(cor_res$p.value, eps = 0.001)

# Create a label for the plot
stats_label <- paste0("Pearson's r = ", r_value, "\n", 
                      "p-value ", ifelse(p_value == "<0.001", "< 0.001", paste("=", p_value)))

# --- 2. GENERATE PLOT ---
p_corr <- ggplot(CARIACO, aes(x = Chlorophyll, y = PrimaryProductivity)) +
  # Add points with some transparency to handle overlapping data
  geom_point(alpha = 0.6, color = "darkblue", size = 2) +
  
  # Add a linear regression line with confidence interval
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  
  # Add the correlation stats to the plot area
  # Adjust x and y coordinates based on your data range (here I put it in the top-left)
  annotate("text", 
           x = min(CARIACO$Chlorophyll, na.rm=T), 
           y = max(CARIACO$PrimaryProductivity, na.rm=T), 
           label = stats_label, 
           hjust = 0, vjust = 1, size = 5, fontface = "italic") +
  
  theme_cowplot() +
  labs(x = expression(Chlorophyll~a~(mg~m^{-3})), 
       y = expression(Primary~Productivity~(mg~C~m^{-2}~d^{-1})),
       title = "Chlorophyll vs. Primary Productivity") + scale_x_log10() + scale_y_log10()

print(p_corr)

