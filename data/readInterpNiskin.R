library(tidyverse)
library(oce)

# import interpolation functions
source('data/interpolateData.r')

niskin_ds <- read.csv("data/BCO-DMO/niskin.csv", na.strings="nd")

niskin_ds$date <- paste(niskin_ds$Year,'-',niskin_ds$Month,'-',niskin_ds$Day, sep='')
niskin_ds$date <- as.Date(niskin_ds$date, format="%Y-%m-%d")

# rename depth variable for interpolation
niskin_ds$depth <- niskin_ds$Depth_real


# CLEAN UP DATA
# remove flagged measurements:
get_q_flags <- niskin_ds %>% select(contains("q_")) %>% names()

# convert q_flags to numeric
niskin_ds_2 <- niskin_ds %>% mutate_at(get_q_flags, as.numeric)

# loop over q_flags and remove flagged values:
for(q_flag in get_q_flags){
  variable = sub('..','',q_flag)
  
  flagged = which(niskin_ds_2[[q_flag]]>0, arr.ind=T)
  
  if (variable=="NO3_NO3_USF"){
    #typo in original dataset needs to be fixed here
    variable = "NO3_NO2_USF"
  }
  
  if (sum(flagged,na.rm=T)>0){
    niskin_ds_2[flagged,][[variable]] <- NA}
}
# note the warning on introduced NAs, this is due to 
# the non-standard NA signifier in the original data set


# add merged nutrients as datapoints
rowwiser <- function(a,b){
  return(rowMeans(cbind(a,b), na.rm=TRUE))
}

niskin_ds_3 <- niskin_ds_2 %>%
  # check if NO3_NO2 - NO2 is negative and return NA if so:
  mutate(NO3_USF = replace(NO3_NO2_USF- NO2_USF, which(NO3_NO2_USF- NO2_USF<0), NA)) %>%
  mutate(NO3_merged = rowwiser(NO3_UDO, NO3_USF),
         PO4_merged = rowwiser(PO4_UDO, PO4_USF),
         SiO4_merged = rowwiser(SiO4_UDO, SiO4_USF))

# create list of numerical variables in datafram to interpolate over
niskin_numeric = c( 
  'O2_umol_kg',
  'O2_ml_L',
  'NO3_UDO',
  'PO4_UDO',
  'SiO4_UDO',
  'NH4_USF',
  'NO2_USF',
  'NO3_NO2_USF',
  'NO3_USF',
  'PO4_USF',
  'SiO4_USF',
  'NO3_merged',
  'PO4_merged',
  'SiO4_merged',
  'pH_corrected',
  'Salinity_bottles',
  'Temperature',
  'Sigma_t',
  'PrimaryProductivity',
  'Chlorophyll',
  'Phaeopigments')

# Export raw data values before interpolation
niskin_raw_out <- niskin_ds_3 %>% select(date, "depth"=Depth_target, all_of(niskin_numeric))
saveRDS(niskin_raw_out, "data/processed/Niskin_RAW.RDS")

# INTERPOLATE
getNiskinIntDepth <- function(depth_from=0, depth_to=100, noofNA=40){
  niskin_temp_store = list()
  
  for (variable in niskin_numeric) {
    # interpolation algorithm: oce-rr
    niskin_temp_store[[variable]] <- interpolateData(niskin_ds_3, variable, depth_from=depth_from, depth_to=depth_to, noofNA = noofNA)
    names(niskin_temp_store[[variable]])[1] <- variable
  }
  
  niskin_ds_cleaned <- niskin_temp_store %>% 
    reduce(full_join, by = "date")
  
  return(niskin_ds_cleaned)
}

niskin_ds_cleaned = getNiskinIntDepth()

# Export and save data
saveRDS(niskin_ds_cleaned, "data/processed/Niskin_qchecked_100m.RDS")
