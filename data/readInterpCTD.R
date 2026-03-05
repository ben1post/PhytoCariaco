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
ctd_temp_int <- interpolateDF(prepdataframe(ctd_ds, "temp"))

# Extract 21 Degree Isotherm
iso21_depth <- ctd_temp_int %>%
  group_by(date) %>%
  filter(depth > 6) %>%
  mutate(iso21 = value_int < 21) %>%
  filter(iso21) %>%
  slice(1)

# Clean up dataframe of Isotherm
iso21_df <- iso21_depth %>%
  rename(Isotherm_21 = depth) %>%
  select(date, Isotherm_21)

# MLD from sigma_t
ctd_sigma_t_int <- interpolateDF(prepdataframe(ctd_ds, "sigma_t"))

ctd_sigma_t_diff <- ctd_sigma_t_int %>%
  group_by(date) %>%
  do(data.frame(sigma_t = .$value_int, sigma_t_diff = c(NA, diff(.$value_int)), depth = .$depth))

# Criterion 2 - relative difference from surface value
mld_depth_2 <- ctd_sigma_t_diff %>%
  group_by(date) %>%
  mutate(sigma_t_surface = sigma_t[depth == min(depth)]) %>%
  filter(depth > 9) %>%
  mutate(mld = abs(sigma_t - sigma_t_surface) >= 0.2) %>%
  filter(mld == TRUE) %>%
  slice(1)

# Clean up MLD data frame
mld_df <- mld_depth_2 %>%
  rename(MLD = depth) %>%
  select(date, MLD)

# Calculate SST from Temperature at surface
sst_10m <- ctd_temp_int %>%
  group_by(date) %>%
  filter(depth <= 10) %>%
  summarize(sst_10m = mean(value_int), .groups = "drop")

# Combine Isotherm and MLD data frames
CTD_combined_data <- list(iso21_df, mld_df, sst_10m) %>% 
  reduce(left_join, by = "date") %>% 
  as.data.frame()

# Add time_month column for merging
CTD_combined_data$time_month <- format(CTD_combined_data$date, format = "%m-%Y")

# Average duplicate months (e.g., 2012-11 had two measurements)
ctd_ds_out <- CTD_combined_data %>% 
  group_by(time_month) %>% 
  summarize(
    Isotherm_21 = mean(Isotherm_21, na.rm = TRUE), 
    MLD = mean(MLD, na.rm = TRUE), 
    sst_10m = mean(sst_10m, na.rm = TRUE),
    .groups = "drop"
  )

# Export data
saveRDS(ctd_ds_out, "data/processed/CTD_Isotherm21_MLD.rds")
