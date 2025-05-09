library(tidyverse)
library(marmap)
library(rnaturalearth)
library(sf)
library(mapsf)

library(ggrepel)
library(RColorBrewer)
library(ggspatial)
# script requires devtools (install.packages("devtools")) to install rnaturalearthhires data for marmap



# coordinates for map bounds
blat1 = -67
blat2 = -63.5
blon1 = 9.5
blon2 = 12


# Get bathymetric data
bat <- getNOAA.bathy(blat1, blat2, blon1, blon2, res = 1, keep = TRUE, path="plots/")
bat_xyz <- as.xyz(bat)
bathy <- as_tibble(bat_xyz)
bathy2 <- bathy %>% mutate(depth_bins = cut(V3, breaks = c(Inf, 0, -50, -100, -200, -500, -1000, 
                                                           -1500, -2000, -Inf))) 

# Import country data
country <- ne_countries(country="Venezuela", scale = "large", returnclass = "sf")

citiez <- ne_download(scale = 10, type = 'populated_places', category = 'cultural') 
citiez_cropped <- st_crop(st_as_sf(citiez), xmin = -69, xmax = -60,
                          ymin = 5, ymax = 15)


# Plot using ggplot and sf
MAP_car <- ggplot() + 
  geom_raster(data = bathy2, aes(x = V1, y = V2, fill = depth_bins), interpolate = TRUE) +
  
  
  scale_fill_manual(values = rev(c("white", brewer.pal(8, "Blues"))), name='Water Depth [m]', 
                    labels=c(">2000","2000", "1500","1000","500","200","100","50","0",""), guide=guide_legend(reverse=TRUE)) +   #, guide = "none") +  
  guides(colour = guide_legend(reverse=TRUE))+
  
  
  # add to map data
  geom_sf(data = country, fill= "antiquewhite", color = "darkblue") +
  
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  
  geom_sf(data = citiez_cropped, color='grey') +
  geom_text_repel(data=citiez_cropped, aes(x=LONGITUDE,y=LATITUDE, label=NAME), hjust = 0, nudge_x = 0.05)+
  
  geom_text(aes(x=-64.5,y=9.7,label="Venezuela"), size=7, color='grey') +
  
  # add mooring:
  geom_point(aes(x=-64.57 ,y=10.5), color='red', size=2) +
  geom_text_repel(aes(x=-64.57 ,y=10.5, label='CARIACO Mooring'), color='red', size=5, hjust = 0.7, nudge_y = 0.2, min.segment.length = 0,
                  arrow = arrow(length = unit(0.015, "npc")), point.padding=1) +
  
  coord_sf(xlim = c(blat1, blat2), 
           ylim = c(blon2, blon1), expand=FALSE, label_axes = list(bottom = "E", top="E", left="N", right="N")) +
  labs(x = "Longitude", y = "Latitude", fill = "Depth (m)") +
  theme_cowplot()#+

#scale_x_continuous(sec.axis = dup_axis()) +
#scale_y_continuous(sec.axis = dup_axis())



options(repr.plot.width=10, repr.plot.height=10)
MAP_car
#ggsave("Map_CAR.pdf", MAP_car, width=10, height=10)



# get globe map with dot for cariaco location

#options(repr.plot.width=3, repr.plot.height=3)
#mtq <- mf_get_mtq()

#mf_export(mtq, filename = "tinyworld.svg", width=3, height=3)

mf_worldmap(lon=-64.57, lat=10.5, pch = 16, lwd = 3, cex = 1,
            water_col = "lightblue", land_col = "antiquewhite",
            border_col = "grey", border_lwd = 1)

#dev.off()