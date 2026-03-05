library(tidyverse)
library(viridis)
library(egg)
library(cowplot)

# --- 1. LOAD DATA ---
# Assuming files exist in your path
NiskinDepth_monthly <- readRDS("data/processed/NiskinDepthIntervals.RDS") 

# --- 2. DEFINE SEASONS & DEPTHS ---
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

# Apply transformations
NiskinDepth_monthly <- NiskinDepth_monthly %>%
  mutate(
    season = get_season_3cat(date),
    # Ensure depth factor order is correct
    depth_rev = factor(depth, levels=c("0-25","25-50","50-75","75-100"))
  )

# --- 3. PLOTTING ---

# Shared Theme
my_theme <- theme_cowplot() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none" # Default to no legend
  )



# Plot 2: Nitrate (NO3_merged)
p_no3 <- ggplot(data = NiskinDepth_monthly, aes(x = season, y = NO3_merged, fill = depth_rev)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8) +
  scale_fill_viridis_d(direction = -1, name = "Depth") + 
  ylab(expression(NO[3]~(mu*M))) +
  my_theme

# Plot 3: Chlorophyll 
# Keeping legend here to serve the whole stack
p_chl <- ggplot(data = NiskinDepth_monthly, aes(x = season, y = Chlorophyll, fill = depth_rev)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8) +
  scale_fill_viridis_d(direction = -1, name = "Depth") + 
  scale_y_log10() + # Log scale often useful for Chl
  ylab(expression(Chl~a~(mg~m^{-3}))) +
  my_theme +
  theme(legend.position = "right") # Enable legend here

# --- 4. ARRANGE PANELS ---
# egg::ggarrange keeps panels aligned even if one has a legend
final_plot <- egg::ggarrange(p_no3, p_chl, ncol = 1)


##### PLOT TWO - Depth X Season X Cluster #####

library(tidyverse)
library(viridis)
library(egg)
library(cowplot)

# --- 1. LOAD DATA ---
NiskinDepth_monthly <- readRDS("data/processed/NiskinDepthIntervals.RDS")

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
NiskinDepth_monthly <- NiskinDepth_monthly %>%
  mutate(
    date = as.Date(date),
    depth_rev = factor(depth, levels = c("0-25", "25-50", "50-75", "75-100")),
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

# --- 4. PLOTTING ---
my_theme <- theme_cowplot() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )

# Plot 1: Nitrate
p_no3 <- ggplot(NiskinDepth_monthly, aes(x = depth_rev, y = NO3_merged, fill = season)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Season") +
  facet_wrap(~ cluster) +
  ylab(expression(NO[3]~(mu*M))) +
  my_theme

# Plot 2: Chlorophyll
p_chl <- ggplot(NiskinDepth_monthly, aes(x = depth_rev, y = Chlorophyll, fill = season)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Season") +
  facet_wrap(~ cluster) +
  scale_y_log10() +
  ylab(expression(Chl~a~(mg~m^{-3}))) +
  my_theme +
  theme(legend.position = "right")

# --- 5. ARRANGE PANELS ---
final_plot <- egg::ggarrange(p_no3, p_chl, ncol = 1)
print(final_plot)




##### PLOT THREE - Season X Cluster ######


library(tidyverse)
library(cowplot)
library(egg)

# --- 1. LOAD & PREP DATA ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

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

# --- 4. DEFINE COLORS & THEME ---
cluster_cols <- c("Cluster 1" = "blue", "Cluster 2" = "red")

my_theme <- theme_cowplot() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )

# --- 5. PLOTTING ---

# Panel 1: Isotherm 21 (inverted y-axis)
p_iso <- ggplot(CARIACO, aes(x = season, y = Isotherm_21, fill = cluster)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = "dodge") +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  scale_y_reverse() +
  ylab("Iso 21 Depth (m)") +
  my_theme

# Panel 2: Nitrate
p_no3 <- ggplot(CARIACO, aes(x = season, y = NO3_merged, fill = cluster)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = "dodge") +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  ylab(expression(NO[3]~(mu*M))) +
  my_theme

# Panel 3: Chlorophyll
p_chl <- ggplot(CARIACO, aes(x = season, y = Chlorophyll, fill = cluster)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = "dodge") +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  ylab(expression(Chl~a~(mg~m^{-3}))) +
  my_theme +
  theme(legend.position = "right") + scale_y_log10()

# --- 6. ARRANGE PANELS ---
final_plot <- egg::ggarrange(p_iso, p_no3, p_chl, ncol = 1)
print(final_plot)
