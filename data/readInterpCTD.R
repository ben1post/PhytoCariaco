library(tidyverse)
library(oce)

# import interpolation functions
source('data/interpolateData.r')

# Read data
ctd_ds <- read.csv("data/BCO-DMO/ctd.csv", na.strings="nd")

# Format date column
ctd_ds$date <- paste(ctd_ds$Year,'-',ctd_ds$Month,'-',ctd_ds$Day, sep='')
ctd_ds$date <- as.Date(ctd_ds$date, format="%Y-%m-%d")

# Interpolate temperature from CTD
ctd_temp_int = interpolateDF(prepdataframe(ctd_ds, "temp"))

# Extract 21 Degree Isotherm
iso21_depth <- ctd_temp_int %>%
  group_by(date) %>%
  filter(depth > 6) %>%
  mutate(iso21 = value_int < 21) %>% # create new column that gives "True" for values at MLD
  filter(iso21 == T) %>% # only take "True" values 
  slice(1) # takes the first occurrence

# Clean up dataframe of Isotherm
iso21_df <- iso21_depth %>%
  rename(Isotherm_21 = depth) %>%
  select(date, Isotherm_21)

# Extract and calculate MLD
# MLD from sigma_t
ctd_sigma_t_int = interpolateDF(prepdataframe(ctd_ds, "sigma_t"))

ctd_sigma_t_diff <- ctd_sigma_t_int %>%
  group_by(date) %>%
  do(data.frame( sigma_t = .$value_int, sigma_t_diff = c(NA,diff(.$value_int)), depth = .$depth))

# Criterion 1 - absolute difference/change in sigma t
mld_depth <- ctd_sigma_t_diff %>%
  group_by(date) %>%
  filter(depth > 9) %>%
  mutate(mld = sigma_t_diff >= 0.125 | sigma_t_diff <= -0.125) %>% # create new column that gives "True" for values at MLD
  filter(mld == T) %>% # only take "True" values 
  slice(1) # takes the first occurrence

# Criterion 2 - relative difference from surface value (This one is used!)
mld_depth_2 <- ctd_sigma_t_diff %>%
  group_by(date) %>%
  filter(depth > 9) %>%
  mutate(mld = sigma_t >= sigma_t[1]+0.2 | sigma_t <= sigma_t[1]-0.2) %>% # create new column that gives "True" for values at MLD
  filter(mld == T) %>% # only take "True" values 
  slice(1) # takes the first occurrence

# Clean up MLD data frame
mld_df <- mld_depth_2 %>%
  rename(MLD = depth) %>%
  select(date, MLD)

# Combine Isotherm and MLD data frames
CTD_combined_data <- list(iso21_df, mld_df) %>% 
  reduce(left_join, by = "date") %>% as.data.frame()


CTD_combined_data$time_month = format(CTD_combined_data$date, format="%m-%Y")

# in 2012-11 there were two measurements, so I need to average these two:
ctd_ds <- CTD_combined_data %>% group_by(time_month) %>% 
  summarize(Isotherm_21 = mean(Isotherm_21), MLD= mean(MLD)) %>%
  mutate(Isotherm_21_lag1=lag(Isotherm_21), Isotherm_21_lag2=lag(Isotherm_21, n=2), Isotherm_21_lag3=lag(Isotherm_21, n=3), Isotherm_21_lag4=lag(Isotherm_21, n=4), Isotherm_21_lag5=lag(Isotherm_21, n=5), Isotherm_21_lag6=lag(Isotherm_21, n=6)) %>%
  ungroup()

# Export data
saveRDS(ctd_ds, "data/processed/CTD_Isotherm21_MLD.rds")
