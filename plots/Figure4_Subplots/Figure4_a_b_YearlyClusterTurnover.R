library(tidyverse)
library(viridis)

require(cowplot)
library(ggdendro)
library(dendextend)

library(cowplot)
library(ggpubr)
library(egg)

library(codyn)

# read interpolated phytoplankton cell count data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

unique(phyto_counts$FuncGroup)


# YEARLY CLUSTERING
ds_genus_yearly <- phyto_counts %>% 
  mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%  
  group_by(Genus, year) %>%
  summarise(Total = sum(counts, na.rm=T))  %>%
  arrange(year)
tail(ds_genus_yearly)

ds_genus_yearly = ds_genus_yearly %>% filter(year>=1996 & year<=2016)

Mesh_genus_yearly <- pivot_wider(ds_genus_yearly, names_from = Genus, values_from = Total, values_fn = function(x) sum(x, na.rm = TRUE), values_fill = 0.0)

Jaccard2_gen <- vegdist(Mesh_genus_yearly[-1], method="jaccard", binary=T)
JM2_gen <- as.matrix(Jaccard2_gen)
colnames(JM2_gen) = Mesh_genus_yearly$year
rownames(JM2_gen) = Mesh_genus_yearly$year


x = hclust(dist(JM2_gen))#, method="ward.D2")

dend <- as.dendrogram(x)

k_3 <- cutree(dend,k = 2, order_clusters_as_data = FALSE) 
nm2 <- setNames( c("Red", "Blue"), c("1","2"))

# export cluster year data for use in coloring other plots:
saveRDS(k_3, "plots/Figure4_Subplots/k_3.RDS")

cols = nm2[as.character(k_3[Mesh_genus_yearly$year])]
cols

options(repr.plot.width=5, repr.plot.height=10)

DENDplot <- dend %>%
  #set("branches_k_color", value = c("blue", "red"), k = 2) %>%
  #set("labels_colors", cols, order_value = TRUE) %>%
  set("labels_cex", 1.) %>% # change color 
  
  hang.dendrogram(hang_height = .1) %>% 
  
  set("leaves_pch", 21) %>% 
  # set("leaves_pch", value = c(21,24), k=2) %>% 
  
  set("leaves_pch", 19) %>% 
  set("leaves_col", cols, order_value = TRUE) %>% 
  set("leaves_bg", cols, order_value = TRUE) %>%
  set("branches_lwd", 0.5) 

DENDplot %>% plot(horiz=T)

DendroPLOT <- ggplot(DENDplot, horiz=FALSE, offset_labels=-0.03) + 
  theme_cowplot(font_size=20) + ylab("Height") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.line.x = element_blank(), 
        strip.clip = "off",
        #plot.margin = margin(0, 25, 0, 0)
  )  +scale_y_continuous(limits=c(0., 1.25)) + scale_x_reverse()

DendroPLOT


# calculate Turnover
# HERE ADD TURNOVER
ds_genus_yearly$year <- as.integer(ds_genus_yearly$year)

# Calculate relative total turnover within replicates
total.res <- turnover(df=ds_genus_yearly,  
                      time.var = "year", 
                      species.var = "Genus",
                      abundance.var = "Total")

total.res$group = k_3[as.character(total.res$year)]
total.res$col = nm2[total.res$group]
total.res$year

total.res$date <- as.Date(as.character(total.res$year), format="%Y")-135

total.res$col

datescale <- scale_x_date(date_labels = "%Y", expand=c(0,0), date_minor_breaks = "1 year", 
                          breaks=c(as.Date("1996/1/1"),as.Date("2000/1/1"),as.Date("2005/1/1"),as.Date("2010/1/1"),
                                   as.Date("2015/1/1"),as.Date("2017/1/1")), 
                          guide = guide_axis(minor.ticks = TRUE))

TurnoverPLOT <- ggplot(data=total.res, aes(x=date,y=total)) + geom_line() + geom_point(aes(col=col), size=3)+ theme_cowplot((font_size=20)) +
  scale_colour_manual(values = c("blue", "red")) + ylab("Species Turnover") + xlab("Date [years]") +labs(colour="Group") + guides(colour="none") + datescale
TurnoverPLOT


### Extract dominant genus per Cluster from yearly aggregate community data
extractGroupsYEAR <- function(year1, year2, year3, year4){
  Group1 <- ds_genus_yearly %>%  
    filter(year >= year1 & year <= year2)
  Group2 <- ds_genus_yearly %>%  
    filter(year >= year3 & year <= year4)
  
  Group = rbind(Group1, Group2)
  
  return(Group)
}

GroupRedGENUS = extractGroupsYEAR(year1=1996, year2=2003, year3=2014, year4=2016) #2016-12-31
GroupRedGENUS_1 = extractGroupsYEAR(year1=1996, year2=2003, year3=2018, year4=2018) #2016-12-31
GroupRedGENUS_2 = extractGroupsYEAR(year1=1995, year2=1995, year3=2014, year4=2016) #2016-12-31
GroupBlueGENUS = extractGroupsYEAR(year1=2004, year2=2013, year3=2017, year4=2017)

returnTop5 <- function(x) {
  top5 <- x %>% group_by(Genus) %>% summarize(Full_sum = sum(Total)) %>% slice_max(order_by=Full_sum, n=5) 
  return(top5)
}

returnTop5(GroupRedGENUS)
returnTop5(GroupRedGENUS_1)
returnTop5(GroupRedGENUS_2)
returnTop5(GroupBlueGENUS)