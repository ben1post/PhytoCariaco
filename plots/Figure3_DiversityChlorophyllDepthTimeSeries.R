library(viridis)
library(egg)
library(cowplot)

# Read depth interval data for phytoplankton:
PhytoDepth_monthly <- readRDS("data/processed/PhytoplanktonDepthIntervals.RDS") 
PhytoDepth_yearly <- readRDS("data/processed/PhytoplanktonDepthIntervals_yrMean.RDS") 

# Test plot of genus richness:
ggplot(data=PhytoDepth_monthly, aes(x=date, y=Pielou_gen, col=as.factor(depth)))+
  geom_point()+
  geom_line(data=PhytoDepth_yearly, aes(x=date, y=Pielou_gen, col=as.factor(depth)), size=1.2) +
  
  geom_point(data=PhytoDepth_yearly, aes(x=date, y=Pielou_gen, col=as.factor(depth)), size=2.7) +
  #geom_segment(data=Depth_yrMean, aes(x=date_start, xend=date_end, y=GenusRichness, yend=GenusRichness, col=as.factor(depth)), size=1.1, lineend="round")+
  scale_color_viridis_d() + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + ylab("Genus Richness") 


# Read depth interval data for niskin (Chlorophyll a):
NiskinDepth_monthly <- readRDS("data/processed/NiskinDepthIntervals.RDS") 
NiskinDepth_yearly <- readRDS("data/processed/NiskinDepthIntervals_yrMean.RDS") 


# Remove unrealiable means from years with missing data:
NiskinDepth_yearly_truncated <- NiskinDepth_yearly %>% filter(as.numeric(year)>=1996 & as.numeric(year)<=2015)

PhytoDepth_yearly_truncated <- PhytoDepth_yearly %>% filter(as.numeric(year)>=1996 & as.numeric(year)<=2015)


#######################################################################################

# create a common date x axis scale for plotting
datescale <- scale_x_date(date_labels = "%Y", expand=c(0,0), date_minor_breaks = "1 year", 
                          breaks=c(as.Date("1996/1/1"),as.Date("2000/1/1"),as.Date("2005/1/1"),as.Date("2010/1/1"),
                                   as.Date("2015/1/1"),as.Date("2017/1/1")), 
                          guide = guide_axis(minor.ticks = TRUE))


# Phytoplankton Diversity and Chlorophyll a time series plot
gl_alpha = 0.3
gl_pointsize = 0.9

tim_0 <- ggplot(data=NiskinDepth_monthly, aes(x=date, y=Chlorophyll, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+    
  
  geom_line(data=NiskinDepth_yearly_truncated, aes(x=date, y=Chlorophyll, col=depth), size=1.) +
  geom_point(data=NiskinDepth_yearly_truncated, aes(x=date, y=Chlorophyll, col=depth), size=2.3) +
  
  scale_color_viridis_d() + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank())  + ylab("Chlorophyll") + guides(col="none")+
  datescale + xlab("Date [Years]") + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border() + scale_y_log10()

tim_1 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=GenusRichness, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+
  
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=GenusRichness, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=GenusRichness, col=depth), size=2.3) +
  scale_color_viridis_d()  + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + ylab("Genus Richness") + guides(col="none")+
  datescale + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()

tim_2 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=Shannon_gen, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+
  
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=Shannon_gen, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=Shannon_gen, col=depth), size=2.3) + guides(col="none")+
  scale_color_viridis_d()  + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + ylab("Shannon Index") +labs(colour="Depth")+
  datescale+ background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()

tim_3 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=Pielou_gen, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+    
  
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=Pielou_gen, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=Pielou_gen, col=depth), size=2.3) +
  
  scale_color_viridis_d() + theme_cowplot()  + ylab("Pielou Index") + guides(col="none")+
  datescale + xlab("Date [Years]") + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()


egg::ggarrange(tim_0,tim_1,tim_2,tim_3, ncol=1)

#######################################################################################
# Reverse depth factor for boxplot
NiskinDepth_monthly$depth_rev <- factor(NiskinDepth_monthly$depth, levels=c("0-25","25-50","50-75","75-100"))
PhytoDepth_monthly$depth_rev <- factor(PhytoDepth_monthly$depth, levels=c("0-25","25-50","50-75","75-100"))


# Boxplot seasonal comparison plot:
sea_0 <- ggplot(data=NiskinDepth_monthly, aes(x=factor(season, levels=c("Upwelling","Rainy")), y=Chlorophyll, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + guides(fill="none") + ylab("Chlorophyll") + scale_y_log10()
#facet_wrap(vars(depth), nrow=1)

sea_1 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=GenusRichness, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank())+ guides(fill="none")  + ylab("Genus Richness")
#facet_wrap(vars(depth), nrow=1)

sea_2 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Shannon_gen, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot()+ theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + guides(fill="none") + labs(colour="Depth") + ylab("Shannon Index")


sea_3 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Pielou_gen, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot()  + xlab("Season") + ylab("Pielou Index")

egg::ggarrange(sea_0, sea_1,sea_2,sea_3)


#######################################################################################

# EXPORT FINAL PLOT:

#options(repr.plot.width=10, repr.plot.height=8)
pdf("plots/exports/Figure3_DepthDivTimeSeries_v4_revResponse.pdf", width=10, height=8)
egg::ggarrange(tim_0,sea_0,tim_1,sea_1, tim_2,sea_2,tim_3, sea_3, 
               ncol = 2, 
               labels = c("a", "b","c","d","e","f","g","h"),
               widths=c(2,0.5)
)

#ggarrange(sea_0,sea_1,sea_2,sea_3, ncol = 1, labels = c("e","f","g","h"))
dev.off()

















############################# V2 Sec Up ########################################


library(tidyverse)
library(viridis)
library(egg)
library(cowplot)

# ==============================================================================
# 1. READ DATA
# ==============================================================================

# Read depth interval data for phytoplankton:
PhytoDepth_monthly <- readRDS("data/processed/PhytoplanktonDepthIntervals.RDS") 
PhytoDepth_yearly <- readRDS("data/processed/PhytoplanktonDepthIntervals_yrMean.RDS") 

# Read depth interval data for niskin (Chlorophyll a):
NiskinDepth_monthly <- readRDS("data/processed/NiskinDepthIntervals.RDS") 
NiskinDepth_yearly <- readRDS("data/processed/NiskinDepthIntervals_yrMean.RDS") 

# Remove unreliable means from years with missing data:
NiskinDepth_yearly_truncated <- NiskinDepth_yearly %>% 
  filter(as.numeric(year) >= 1996 & as.numeric(year) <= 2015)

PhytoDepth_yearly_truncated <- PhytoDepth_yearly %>% 
  filter(as.numeric(year) >= 1996 & as.numeric(year) <= 2015)


# ==============================================================================
# 2. TIME SERIES PLOTS (Left Column)
# ==============================================================================

# create a common date x axis scale for plotting
datescale <- scale_x_date(date_labels = "%Y", expand=c(0,0), date_minor_breaks = "1 year", 
                          breaks=c(as.Date("1996/1/1"), as.Date("2000/1/1"), as.Date("2005/1/1"), 
                                   as.Date("2010/1/1"), as.Date("2015/1/1"), as.Date("2017/1/1")), 
                          guide = guide_axis(minor.ticks = TRUE))

# Phytoplankton Diversity and Chlorophyll a time series plot settings
gl_alpha = 0.3
gl_pointsize = 0.9

tim_0 <- ggplot(data=NiskinDepth_monthly, aes(x=date, y=Chlorophyll, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+    
  geom_line(data=NiskinDepth_yearly_truncated, aes(x=date, y=Chlorophyll, col=depth), size=1.) +
  geom_point(data=NiskinDepth_yearly_truncated, aes(x=date, y=Chlorophyll, col=depth), size=2.3) +
  scale_color_viridis_d() + theme_cowplot() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  ylab("Chlorophyll") + guides(col="none") +
  datescale + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border() + scale_y_log10()

tim_1 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=GenusRichness, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=GenusRichness, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=GenusRichness, col=depth), size=2.3) +
  scale_color_viridis_d() + theme_cowplot() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  ylab("Genus Richness") + guides(col="none") +
  datescale + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()

tim_2 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=Shannon_gen, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=Shannon_gen, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=Shannon_gen, col=depth), size=2.3) + 
  scale_color_viridis_d() + theme_cowplot() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  ylab("Shannon Index") + guides(col="none") +
  datescale + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()

tim_3 <- ggplot(data=PhytoDepth_monthly, aes(x=date, y=Pielou_gen, col=depth))+
  geom_point(alpha=gl_alpha, size=gl_pointsize)+    
  geom_line(data=PhytoDepth_yearly_truncated, aes(x=date, y=Pielou_gen, col=depth), size=1.) +
  geom_point(data=PhytoDepth_yearly_truncated, aes(x=date, y=Pielou_gen, col=depth), size=2.3) +
  scale_color_viridis_d() + theme_cowplot() + 
  ylab("Pielou Index") + guides(col="none") +
  datescale + xlab("Date [Years]") + background_grid(major = "x", minor = "x", size.minor = 0.5) + panel_border()


# ==============================================================================
# 3. DEFINE SEASONS & PREP DATA FOR BOXPLOTS
# ==============================================================================

get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  season_lookup <- c(
    "12" = "Primary", "01" = "Primary", "02" = "Primary", 
    "03" = "Primary", "04" = "Primary", "05" = "Primary",
    "06" = "Secondary", "07" = "Secondary", "08" = "Secondary",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  factor(season_lookup[months], levels = c("Primary", "Secondary", "Rainy"))
}

# Apply seasons and reverse depth factor for Niskin Data
NiskinDepth_monthly <- NiskinDepth_monthly %>%
  mutate(
    season = get_season_3cat(date),
    depth_rev = factor(depth, levels=c("0-25","25-50","50-75","75-100"))
  )

# Apply seasons and reverse depth factor for Phyto Data
PhytoDepth_monthly <- PhytoDepth_monthly %>%
  mutate(
    season = get_season_3cat(date),
    depth_rev = factor(depth, levels=c("0-25","25-50","50-75","75-100"))
  )


# ==============================================================================
# 4. SEASONAL BOXPLOTS (Right Column)
# ==============================================================================

# Y-axis numbers (text) are kept, but the titles (e.g., "Chlorophyll") are removed.

sea_0 <- ggplot(data=NiskinDepth_monthly, aes(x=season, y=Chlorophyll, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) + 
  guides(fill="none") + scale_y_log10()

sea_1 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=GenusRichness, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) + 
  guides(fill="none")

sea_2 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Shannon_gen, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) + 
  guides(fill="none")

sea_3 <- ggplot(data=PhytoDepth_monthly, aes(x=season, y=Pielou_gen, fill=depth_rev))+
  geom_boxplot()+
  scale_fill_viridis_d(direction=-1) + theme_cowplot() + 
  xlab("Season") + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 15, hjust = 1)) 


# ==============================================================================
# 5. EXPORT FINAL PLOT
# ==============================================================================

pdf("plots/exports/Figure3_DepthDivTimeSeries_v4_revResponse.pdf", width=10, height=8)
egg::ggarrange(tim_0, sea_0, 
               tim_1, sea_1, 
               tim_2, sea_2, 
               tim_3, sea_3, 
               ncol = 2, 
               labels = c("a", "b","c","d","e","f","g","h"),
               widths = c(2, 0.8) # Adjusted to give boxplots slightly more width
)
dev.off()
