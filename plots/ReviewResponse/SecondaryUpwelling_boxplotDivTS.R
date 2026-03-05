library(viridis)
library(egg)
library(cowplot)


# --- 1. LOAD DATA ---
PhytoDepth_monthly <- readRDS("data/processed/PhytoplanktonDepthIntervals.RDS") 
NiskinDepth_monthly <- readRDS("data/processed/NiskinDepthIntervals.RDS") 

# --- 2. DEFINE NEW SEASON MAPPING ---
# Logic: 
# Upwelling: Dec-May
# Secondary Upwelling: Jun-Aug
# Rainy: Sep-Nov

get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  
  # Create a named vector for easy lookup
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  
  # Return the season as a factor with specific order
  factor(season_lookup[months], 
         levels = c("Upwelling", "Secondary Upwelling", "Rainy"))
}

# --- 3. APPLY TO DATA FRAMES ---
PhytoDepth_monthly$season <- get_season_3cat(PhytoDepth_monthly$date)
NiskinDepth_monthly$season <- get_season_3cat(NiskinDepth_monthly$date)

# Also ensure depth factors are consistent (as in your original code)
PhytoDepth_monthly$depth_rev <- factor(PhytoDepth_monthly$depth, levels=c("0-25","25-50","50-75","75-100"))
NiskinDepth_monthly$depth_rev <- factor(NiskinDepth_monthly$depth, levels=c("0-25","25-50","50-75","75-100"))


library(viridis)
library(egg)
library(cowplot)

# Define a shared theme element to rotate x-axis labels slightly 
# (helps fit "Secondary Upwelling")
my_theme <- theme_cowplot() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1), # Angle text to fit
    legend.position = "none" # Hide legends for individual plots (cleaner)
  )

# 1. Chlorophyll
sea_0 <- ggplot(data=NiskinDepth_monthly, aes(x=season, y=Chlorophyll, fill=depth_rev))+
  geom_boxplot() +
  scale_fill_viridis_d(direction=-1, name="Depth") + 
  scale_y_log10() +
  ylab("Chlorophyll") +
  my_theme

# 2. Genus Richness
sea_1 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=GenusRichness, fill=depth_rev))+
  geom_boxplot() +
  scale_fill_viridis_d(direction=-1) + 
  ylab("Genus Richness") +
  my_theme

# 3. Shannon Index
sea_2 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Shannon_gen, fill=depth_rev))+
  geom_boxplot() +
  scale_fill_viridis_d(direction=-1) + 
  ylab("Shannon Index") +
  my_theme

# 4. Pielou Index (Keep legend here if you want one legend for the whole figure)
sea_3 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Pielou_gen, fill=depth_rev))+
  geom_boxplot() +
  scale_fill_viridis_d(direction=-1, name="Depth") + 
  ylab("Pielou Index") +
  my_theme +
  theme(legend.position = "right") # Keep legend only on the last one?

# Arrange
# Note: Since sea_3 has a legend and others don't, alignment might be tricky.
# It is often better to extract the legend separately or use 'cowplot::plot_grid'
final_plot <- egg::ggarrange(sea_0, sea_1, sea_2, sea_3, nrow = 2)
