# Script to merge all processed datasets into cohesive data frame to use for the analysis

library(tidyverse)
library(vegan)
# set strings as factors to false
options(stringsAsFactors = FALSE)

# READ ERA DATA
wind_ds <- readRDS("data/processed/ERA5_data.rds")
wind_ds$date = as.Date(wind_ds$time)
wind_ds$time_month = format(wind_ds$date, format="%m-%Y")
wind_ds <- wind_ds %>% select(-date, -time)
str(wind_ds)

wind_ds <- wind_ds %>% mutate(u10_lag1=lag(u10), u10_lag2=lag(u10, n=2), u10_lag3=lag(u10, n=3), u10_lag4=lag(u10, n=4), u10_lag5=lag(u10, n=5), u10_lag6=lag(u10, n=6))



# READ NISKIN DATA
niskin_ds <- readRDS("data/processed/Niskin_qchecked_100m.RDS")
niskin_ds$time_month = format(niskin_ds$date, format="%m-%Y")
niskin_ds <- niskin_ds %>% select(-date) %>%
  mutate(NO3_merged_lag1=lag(NO3_merged), NO3_merged_lag2=lag(NO3_merged, n=2), NO3_merged_lag3=lag(NO3_merged, n=3),
         PO4_merged_lag1=lag(PO4_merged), PO4_merged_lag2=lag(PO4_merged, n=2), PO4_merged_lag3=lag(PO4_merged, n=3),
         SiO4_merged_lag1=lag(SiO4_merged), SiO4_merged_lag2=lag(SiO4_merged, n=2), SiO4_merged_lag3=lag(SiO4_merged, n=3),
         Temperature_lag1=lag(Temperature, n=1), Temperature_lag2=lag(Temperature, n=2), Temperature_lag3=lag(Temperature, n=3),
         Salinity_bottles_lag1=lag(Salinity_bottles, n=1), Salinity_bottles_lag2=lag(Salinity_bottles, n=2), Salinity_bottles_lag3=lag(Salinity_bottles, n=3))
str(niskin_ds)

# READ PHYTOPLANKTON DATA
diversity_ds <- readRDS("data/processed/PhytoplanktonAbundanceDiversity.RDS")
str(diversity_ds)

# READ CTD DATA
ctd_ds <- readRDS("data/processed/CTD_Isotherm21_MLD.rds")
str(ctd_ds)

# READ AMO DATA
amo_table <- read.table("data/AMO/amo_monthly.txt", header=FALSE, skip=1)

amo_df <- amo_table %>% gather(Month, Value, -V1)

nm2 <- setNames( c("01","02","03","04","05","06","07","08","09","10","11","12"), c("V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"))
amo_df$Month = nm2[amo_df$Month]

names(amo_df) <- c("Year","Month","AMO")

amo_ds <- amo_df %>% 
  mutate(time_month = paste(Month, Year, sep="-")) %>%
  select(AMO, time_month) %>%
  
  mutate(AMO_lag1=lag(AMO), AMO_lag2=lag(AMO, n=2), AMO_lag3=lag(AMO, n=3), AMO_lag4=lag(AMO, n=4), AMO_lag5=lag(AMO, n=5), AMO_lag6=lag(AMO, n=6), AMO_lag11=lag(AMO, n=11),AMO_lag12=lag(AMO, n=12))

str(amo_ds)


# READ MEI v.2 DATA
meiv2_lines <- readLines("data/MEIv2/meiv2.data")
meiv2_lines = meiv2_lines[c(-1,-48:-51)] # remove unecessary/text lines

meiv2_table <- read.table(textConnection(meiv2_lines), header=FALSE, stringsAsFactors = FALSE)


meiv2_df <- meiv2_table %>% gather(Month, Value, -V1)
meiv2_df$Month = nm2[meiv2_df$Month]

names(meiv2_df) <- c("Year","Month","MEIv2")

meiv2_ds <- meiv2_df %>% 
  mutate(time_month = paste(Month, Year, sep="-")) %>%
  select(MEIv2, time_month) %>%
  
  mutate(MEIv2_lag1=lag(MEIv2), MEIv2_lag2=lag(MEIv2, n=2), MEIv2_lag3=lag(MEIv2, n=3), MEIv2_lag4=lag(MEIv2, n=4), MEIv2_lag5=lag(MEIv2, n=5), MEIv2_lag6=lag(MEIv2, n=6), MEIv2_lag12=lag(MEIv2, n=12))

str(meiv2_ds)


# COMBINE ALL DATASETS INTO ONE DATA FRAME
CARIACO_dat_joined <- list(wind_ds, 
                           niskin_ds,
                           diversity_ds,
                           ctd_ds,
                           amo_ds,
                           meiv2_ds
) %>% 
  reduce(full_join, by = c("time_month"))

# Remove Years out of scope of CARIACO time series
CARIACO_dat_joined_truncated <- CARIACO_dat_joined[which(CARIACO_dat_joined$time_month=="11-1995"):which(CARIACO_dat_joined$time_month=="01-2017"),]

names(CARIACO_dat_joined_truncated)

# Save and export
saveRDS(CARIACO_dat_joined_truncated, "data/processed/CARIACO_EnvData_combined.rds")
