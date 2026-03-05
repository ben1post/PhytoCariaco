library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(scales)


# read combined dataset:
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")
CARIACO$year <- format(CARIACO$date, format="%Y")

# calculate satistics for all variables in CARIACO dataset:
CARIACOyearly <- CARIACO %>% 
  group_by(year) %>% 
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))) %>% 
  ungroup()

# calculate missing values per year
CARIACOyearly_NA <- CARIACO %>% group_by(year) %>% 
  summarize(across(where(is.numeric), ~sum(is.na(.))))


# Function to get yearly mean
# fill_na: if TRUE, fills missing months with monthly climatology; if FALSE, uses raw data only
getym <- function(variable = '', fill_na = FALSE) {
  
  df <- CARIACO %>% 
    select(date, all_of(variable)) %>% 
    rename(value = all_of(variable))
  
  out <- df %>% 
    mutate(Month = month(date)) %>% 
    group_by(Month) %>%
    mutate(monthly_average = mean(value, na.rm = TRUE))
  
  if (fill_na) {
    # Replace NAs with monthly climatology
    out <- out %>%
      mutate(value_used = replace(value, is.na(value), monthly_average[is.na(value)]))
  } else {
    # Keep NAs as is
    out <- out %>%
      mutate(value_used = value)
  }
  
  out$Year <- as.numeric(format(out$date, format = "%Y"))
  
  df_yearly <- out %>% 
    group_by(Year) %>%
    summarize(
      value_filled = mean(value_used, na.rm = TRUE),  # Uses filled or raw depending on fill_na
      value_raw = mean(value, na.rm = TRUE),
      n_na = sum(is.na(value)),
      .groups = "drop"
    ) %>%
    mutate(value = if (fill_na) value_filled else value_raw) %>%
    select(Year, value, value_raw, n_na)
  
  return(df_yearly)
}

# reformat and scale data variables for visualisation
CARIACO$u10_negative = -CARIACO$u10
#CARIACO$Isotherm_21_negative = -CARIACO$Isotherm_21

CARIACO$Chlorophyll_log = log1p(CARIACO$Chlorophyll)
CARIACO$PrimaryProductivity_log = log1p(CARIACO$PrimaryProductivity)
CARIACO$Abundance_Diatom_log = log1p(CARIACO$Abundance_Diatom)
CARIACO$Abundance_Hapto_log = log1p(CARIACO$Abundance_Hapto)
CARIACO$Abundance_Dino_log = log1p(CARIACO$Abundance_Dino)
CARIACO$Abundance_Cyano_log = log1p(CARIACO$Abundance_Cyano)
CARIACO$Abundance_Nanoflagellate_log = log1p(CARIACO$Abundance_Nanoflagellate)

range01 <- function(x){scale(x)}
#{(x-min(x))/(max(x)-min(x))-0.5}

# Define variables to plot and their display names
vars <- c("MEIv2", "AMO", "u10_negative", "tp", "Isotherm_21", "Salinity_bottles", "sst_10m",
          "NO3_merged", "PO4_merged", "SiO4_merged", "Chlorophyll_log",
          "PrimaryProductivity_log",
          "Abundance_Diatom_log", "Abundance_Hapto_log", "Abundance_Dino_log", 
          "Abundance_Cyano_log", "Abundance_Nanoflagellate_log",
          "GenusRichness", "Shannon_gen", "Pielou_gen")

renamed_vars <- c("MEIv2", "AMO", "Wind speed", "Precipitation", "Isotherm Depth", "Salinity", "SST",
                  "NO3", "PO4", "SiO4", "Chlorophyll a",
                  "Primary Productivity",
                  "Diatoms", "Haptophytes", "Dinoflagellates", "Cyanobacteria", "Nanoflagellates",
                  "Richness", "Shannon", "Pielou")

sources <- c("Climate", "Climate", "Climate", "Climate", "Physical", "Physical", "Physical",
             "Physical", "Physical", "Physical", "Biological",
             "Biological",
             "Biological", "Biological", "Biological", "Biological", "Biological",
             "Diversity", "Diversity", "Diversity")

dat_temp_store = list()

# Set this to TRUE to fill NA months with monthly means across the entire time series
use_filled_na <- FALSE

for (i in 1:length(vars)) {
  dat_temp_store[[vars[i]]] <- getym(vars[i], fill_na = use_filled_na)
  dat_temp_store[[vars[i]]]$zscore_raw <- scale(dat_temp_store[[vars[i]]]$value_raw)
  dat_temp_store[[vars[i]]]$zscore <- scale(dat_temp_store[[vars[i]]]$value)
  dat_temp_store[[vars[i]]]$var <- renamed_vars[i]
  dat_temp_store[[vars[i]]]$source <- sources[i]
}

library(data.table)
full_yearly_dat <- rbindlist(dat_temp_store)
head(full_yearly_dat)



yearly_zscore = full_yearly_dat #%>% filter(year>=1996 & year<=2015)

var_reorder <- c("MEIv2", "AMO", "Wind speed", "Precipitation", "SST", "Isotherm Depth",
                 "Salinity",
                 "NO3", "PO4", "SiO4", 
                 "Chlorophyll a", "Primary Productivity", "Diatoms", "Haptophytes",
                 "Dinoflagellates", "Cyanobacteria",
                 "Nanoflagellates",
                 "Richness", "Shannon", "Pielou")
yearly_zscore$ord.var <- factor(yearly_zscore$var, ordered = TRUE, levels = rev(var_reorder))

str(yearly_zscore)


panel_reorder = c("Climate","Physical","Biological","Diversity")
yearly_zscore$ord.source <- factor(yearly_zscore$source, ordered=TRUE, levels = panel_reorder)
str(yearly_zscore)


cols <- rev(brewer.pal(n = 5, name = "RdBu"))



# --- PREPARE DATA FOR MISSING DATA INDICATORS ---

# Create x-axis position mapping (years 1996-2015)
year_levels <- as.character(1996:2015)

# Compute y positions within each facet panel
y_pos_lookup <- yearly_zscore %>%
  select(var, ord.var, ord.source) %>%
  distinct() %>%
  group_by(ord.source) %>%
  arrange(ord.var) %>%
  mutate(y_pos = row_number()) %>%
  ungroup() %>%
  select(var, ord.source, y_pos)

# Join positions to main data and compute x positions
yearly_zscore_pos <- yearly_zscore %>%
  left_join(y_pos_lookup, by = c("var", "ord.source")) %>%
  mutate(x_pos = match(as.character(Year), year_levels))

# Filter for different missing data categories
df_na_two <- yearly_zscore_pos %>% filter(n_na == 2)
df_na_medium <- yearly_zscore_pos %>% filter(n_na == 3)
df_na_high <- yearly_zscore_pos %>% filter(n_na >= 4)
df_na_border <- yearly_zscore_pos %>% filter(n_na > 2)



# PLOT:

options(repr.plot.width=14, repr.plot.height=10)

ggplot(yearly_zscore, aes(x=as.character(Year), y=ord.var, fill=zscore))+ geom_tile()+
  # geom_tile(data=yearly_zscore %>% filter(n_na==1), aes(x=as.character(Year), y=ord.var, fill=zscore_raw), color="black", size=.3, width=0.97, height=0.97) +
  # geom_tile(data=yearly_zscore %>% filter(n_na==2), aes(x=as.character(Year), y=ord.var, fill=zscore_raw), color="blue", size=.5, width=0.95, height=0.95) +
  # geom_tile(data=yearly_zscore %>% filter(n_na==3), aes(x=as.character(Year), y=ord.var, fill=zscore_raw), color="green", size=.8, width=0.94, height=0.94) +
  # geom_tile(data=yearly_zscore %>% filter(n_na>=4), aes(x=as.character(Year), y=ord.var, fill=zscore_raw), color="red", size=.8, width=0.94, height=0.94) +
  # 
  # Black border for cells with >2 missing
  # Small square for cells with exactly 2 missing
  geom_point(
    data = df_na_two,
    aes(x = as.character(Year), y = ord.var),
    shape = 0,
    size = 0.1,
    stroke = 0.8,
    color = "black"
  ) +
  geom_tile(
    data = df_na_border,
    aes(x = as.character(Year), y = ord.var),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  
  # Single diagonal (>2 and <4 missing)
  geom_segment(
    data = df_na_medium,
    aes(x = x_pos - 0.5, xend = x_pos + 0.5,
        y = y_pos - 0.5, yend = y_pos + 0.5),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  
  # Crossed diagonals (>=4 missing)
  geom_segment(
    data = df_na_high,
    aes(x = x_pos - 0.5, xend = x_pos + 0.5,
        y = y_pos - 0.5, yend = y_pos + 0.5),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    data = df_na_high,
    aes(x = x_pos - 0.5, xend = x_pos + 0.5,
        y = y_pos + 0.5, yend = y_pos - 0.5),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  
  facet_grid(ord.source ~ ., scales = "free_y", space = "free_y") +
  
  #coord_fixed()+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 10)) +  
  #scale_fill_gradient2(low = "#1a2c87",
  #                   mid = "#FFFFFF",
  #                   high = "#b11e23") +
  
  #scale_fill_distiller(palette = "RdBu") +
  scale_fill_gradientn(colours = cols, 
                       values = rescale(c(-3, -0.5, 0, 0.5, 4)),
                       guide = "colorbar", limits=c(-3, 4))+
  
  ggtitle("EnvData - Yearly Avg zscores")+
  labs(x="Year", fill="z-score") + theme_cowplot(font_size=20)+
  theme(
    axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
    axis.title.y=element_blank(),
    plot.title = element_text(size = 15, face = "bold"),
    axis.line=element_blank())+
  xlim(as.character(1996:2015))


#ggsave("plots/exports/Figure2_ZScores_v4_revresponse.pdf", width=14, height=10)



