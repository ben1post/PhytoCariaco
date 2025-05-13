library(tidyverse)

niskin_raw <- readRDS("data/processed/Niskin_RAW.RDS")

# filter unnecessary depths out of niskin
niskin_raw_filtdep <- niskin_raw %>% filter(!depth>100)
dim(niskin_raw_filtdep)
# remove outliers
niskin_raw_filtdep <- niskin_raw %>% filter(!Chlorophyll>20)
dim(niskin_raw_filtdep)
# add month column
niskin_raw_filtdep$month <- format(niskin_raw_filtdep$date, format="%m")

# add season column
#nm2 <- setNames( rep(c("Winter", "Spring", "Summer", "Fall"),
#        each = 3), c("12","01","02","03","04","05","06","07","08","09","10","11"))
nm2 <- setNames( c(rep(c("Upwelling", "Rainy"),
                       each = 5),"Rainy","Rainy"), c("12","01","02","03","04","05","06","07","08","09","10","11"))
#add columns to data frame 
niskin_raw_filtdep$year = format(niskin_raw_filtdep$date, "%Y")
niskin_raw_filtdep$season = nm2[niskin_raw_filtdep$month]


# import interpolation functions
source('data/interpolateData.r')

niskin_vars = c("Chlorophyll","PrimaryProductivity", "NO3_merged")

# interpolation function
getNiskinInterpCounts <- function(depth_from=0, depth_to=100, noofNA=20){
  niskin_temp_store = list()
  
  for (variable in niskin_vars) {
    # interpolation algorithm: oce-rr
    niskin_temp_store[[variable]] <- interpolateData(niskin_raw_filtdep, variable, depth_from=depth_from, depth_to=depth_to, noofNA=noofNA, int_func='unesco')
    names(niskin_temp_store[[variable]])[1] <- variable
  }
  
  niskin_ds_cleaned <- niskin_temp_store %>% 
    reduce(full_join, by = "date") %>%
    mutate("depth"=paste(depth_from, depth_to, sep="-"))
  
  
  return(niskin_ds_cleaned)
}

# calculate depth discrete intervals
niskin_int_0 <- getNiskinInterpCounts(depth_from=0, depth_to=25, noofNA=5)
niskin_int_1 <- getNiskinInterpCounts(depth_from=25, depth_to=50, noofNA=5)
niskin_int_2 <- getNiskinInterpCounts(depth_from=50, depth_to=75, noofNA=5)
niskin_int_3 <- getNiskinInterpCounts(depth_from=75, depth_to=100, noofNA=5)
# combine into single data frame
niskin_dat <- rbind(niskin_int_0, niskin_int_1, niskin_int_2, niskin_int_3)
niskin_dat$depth <- factor(niskin_dat$depth, levels=c("75-100","50-75","25-50", "0-25"))

# add month column
niskin_dat$month <- format(niskin_dat$date, format="%m")

# add season column
#nm2 <- setNames( rep(c("Winter", "Spring", "Summer", "Fall"),
#        each = 3), c("12","01","02","03","04","05","06","07","08","09","10","11"))
nm2 <- setNames( c(rep(c("Upwelling", "Rainy"),
                       each = 6)), c("12","01","02","03","04","05","06","07","08","09","10","11"))
#add year and season columns to data frame 
niskin_dat$year = format(niskin_dat$date, "%Y")
niskin_dat$season = nm2[niskin_dat$month]

names(niskin_dat)
# calculate yearly means
niskin_yrMean <- niskin_dat %>% 
  group_by(depth,year) %>%
  summarize(date_start=min(date), date_end=max(date), date=mean(date), 
            Chlorophyll=mean(Chlorophyll,na.rm=T), PrimaryProductivity=mean(PrimaryProductivity,na.rm=T))


# Save depth discrete data
saveRDS(niskin_dat, "data/processed/NiskinDepthIntervals.RDS")
saveRDS(niskin_yrMean, "data/processed/NiskinDepthIntervals_yrMean.RDS")
