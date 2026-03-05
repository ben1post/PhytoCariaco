library(tidyverse)
library(cowplot)
library(egg)

# Load data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

CARIACO$u10_negative = -CARIACO$u10

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

# Prep data
CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    season = get_season_3cat(date),
    AMO_phase = factor(ifelse(AMO >= 0, "AMO+", "AMO-"), levels = c("AMO-", "AMO+"))
  )

# Shared theme
my_theme <- theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        legend.position = "none")

amo_cols <- c("AMO-" = "steelblue", "AMO+" = "tomato")

# Individual plots
p_wind <- ggplot(CARIACO, aes(x = season, y = u10_negative, fill = AMO_phase)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = amo_cols) +
  ylab("Wind Speed") + my_theme

p_iso <- ggplot(CARIACO, aes(x = season, y = Isotherm_21, fill = AMO_phase)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = amo_cols) +
  ylab("Isotherm Depth (m)") + my_theme + scale_y_reverse()

p_no3 <- ggplot(CARIACO, aes(x = season, y = NO3_merged, fill = AMO_phase)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = amo_cols) +
  ylab(expression(NO[3]~(mu*M))) + my_theme

p_chl <- ggplot(CARIACO, aes(x = season, y = Chlorophyll, fill = AMO_phase)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = amo_cols) +
  scale_y_log10() +
  ylab(expression(Chl~a~(mg~m^{-3}))) + my_theme

p_diatom <- ggplot(CARIACO, aes(x = season, y = Abundance_Diatom, fill = AMO_phase)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = amo_cols) +
  scale_y_log10() +
  ylab("Diatom Abundance") + my_theme +
  theme(legend.position = "right") +
  labs(fill = "AMO Phase")

# Combine
egg::ggarrange(p_wind, p_iso, p_no3, p_chl, p_diatom, ncol = 5)



####### double check cluster effects #######


library(tidyverse)
library(cowplot)
library(egg)


# Define seasons
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

CARIACO <- CARIACO %>%
  mutate(
    date = as.Date(paste(time_month, "-15", sep = ""), format = "%m-%Y-%d"),
    season = get_season_3cat(date),
    cluster = case_when(
      (date >= as.Date("1996-01-01") & date <= as.Date("2004-12-31")) |
        (date >= as.Date("2014-01-01") & date <= as.Date("2017-12-31")) ~ "Cluster 1",
      (date >= as.Date("2005-01-01") & date <= as.Date("2013-12-31")) ~ "Cluster 2",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cluster))

# Shared theme
my_theme <- theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        legend.position = "none")

cluster_cols <- c("Cluster 1" = "blue", "Cluster 2" = "red")

# Individual plots
p_wind <- ggplot(CARIACO, aes(x = season, y = u10_negative, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = cluster_cols) +
  ylab("Wind Speed") + my_theme

p_iso <- ggplot(CARIACO, aes(x = season, y = Isotherm_21, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_reverse() +
  ylab("Isotherm Depth (m)") + my_theme

p_no3 <- ggplot(CARIACO, aes(x = season, y = NO3_merged, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = cluster_cols) +
  ylab(expression(NO[3]~(mu*M))) + my_theme

p_chl <- ggplot(CARIACO, aes(x = season, y = Chlorophyll, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_log10() +
  ylab(expression(Chl~a~(mg~m^{-3}))) + my_theme

p_diatom <- ggplot(CARIACO, aes(x = season, y = Abundance_Diatom, fill = cluster)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = cluster_cols) +
  scale_y_log10() +
  ylab("Diatom Abundance") + my_theme +
  theme(legend.position = "right") +
  labs(fill = "Cluster")

# Combine
egg::ggarrange(p_wind, p_iso, p_no3, p_chl, p_diatom, ncol = 5)

