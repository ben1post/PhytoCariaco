library(tidyverse)
library(vegan)

Mesh_phyInt_genus<-readRDS("data/processed/PhytoplanktonRawGenusCounts.RDS")

# import interpolation functions
source('data/interpolateData.r')

# get list of Genera to iterate over:
phyto_Genus = names(Mesh_phyInt_genus[,c(-1,-2)])

# Function to interpolate raw counts and calculate diversity metrics for depth intervals
getPhytoInterpCounts <- function(depth_from=0, depth_to=100, noofNA=20){
  phyto_temp_store = list()
  
  for (variable in phyto_Genus) {
    # interpolation algorithm: oce-rr
    phyto_temp_store[[variable]] <- interpolateData(Mesh_phyInt_genus, variable, depth_from=depth_from, depth_to=depth_to, noofNA=noofNA, int_func='unesco')
    names(phyto_temp_store[[variable]])[1] <- variable
  }
  
  phyto_ds_cleaned <- phyto_temp_store %>% 
    reduce(full_join, by = "date") %>% na.omit()
  
  phyto_pivot <- phyto_ds_cleaned %>% 
    pivot_longer(cols=-date, names_to = "Genus", values_to = "counts") 
  
  phyto_counts = phyto_pivot[complete.cases(phyto_pivot$counts),]
  
  # get Genus Matrix:
  ds_genus <- phyto_counts %>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts))  %>%
    arrange(date)    
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  Mesh_genus = Mesh_genus[rowSums(Mesh_genus[, -1])>0, ]
  
  Shannon_gen <- diversity(Mesh_genus[-1])
  Pielou_gen <- Shannon_gen/log(specnumber(Mesh_genus[-1]))
  GenusRichness <- apply(Mesh_genus[-1]>0,1,sum)
  
  diversity_ds = data.frame("date"=Mesh_genus[[1]], "depth"=paste(depth_from, depth_to, sep="-"), GenusRichness, Shannon_gen, Pielou_gen)
  return(diversity_ds)    
}

# interpolate across 4 depth intervals:
phyto_int_0 <- getPhytoInterpCounts(depth_from=0, depth_to=25, noofNA=5)
phyto_int_1 <- getPhytoInterpCounts(depth_from=25, depth_to=50, noofNA=5)
phyto_int_2 <- getPhytoInterpCounts(depth_from=50, depth_to=75, noofNA=5)
phyto_int_3 <- getPhytoInterpCounts(depth_from=75, depth_to=100, noofNA=5)

# combine into a single data frame
Depth_Dat <- rbind(phyto_int_0,phyto_int_1,phyto_int_2,phyto_int_3)
# convert depth to factor with correct level order
Depth_Dat$depth <- factor(Depth_Dat$depth, levels=c("75-100","50-75","25-50", "0-25"))

names(Depth_Dat)
# convert date
Depth_Dat$date <- as.Date(Depth_Dat$date)

ggplot(data=Depth_Dat, aes(x=date, y=Shannon_gen, col=depth))+
  geom_point()+
  geom_smooth(se=FALSE)+scale_color_viridis_d()

# add month
Depth_Dat$month <- format(Depth_Dat$date, format="%m")

# add season column
nm2 <- setNames( c(rep(c("Upwelling", "Rainy"),
                       each = 6)), c("12","01","02","03","04","05","06","07","08","09","10","11"))
# add columns to data frame 
Depth_Dat$year = format(Depth_Dat$date, "%Y")
Depth_Dat$season = factor(nm2[Depth_Dat$month], levels = c("Upwelling", "Rainy"))

head(Depth_Dat)

# calculate yearly means
Depth_yrMean <- Depth_Dat %>% 
  group_by(depth,year) %>%
  summarize(date_start=min(date), date_end=max(date),date=mean(date), 
            GenusRichness=mean(GenusRichness,na.rm=T), Shannon_gen=mean(Shannon_gen,na.rm=T), Pielou_gen=mean(Pielou_gen,na.rm=T))
# convert depth to factor with correct level order
Depth_yrMean$depth <- factor(Depth_yrMean$depth, levels=c("75-100","50-75","25-50", "0-25"))


# Save depth discrete data
saveRDS(Depth_Dat, "data/processed/PhytoplanktonDepthIntervals.RDS")
saveRDS(Depth_yrMean, "data/processed/PhytoplanktonDepthIntervals_yrMean.RDS")
