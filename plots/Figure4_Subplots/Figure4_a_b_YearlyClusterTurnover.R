library(tidyverse)
library(viridis)

require(cowplot)
library(ggdendro)
library(dendextend)

library(cowplot)
library(ggpubr)
library(egg)

library(vegan)
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

x = hclust(as.dist(JM2_gen))

dend <- as.dendrogram(x)

# Updated to 3 clusters
k_3 <- cutree(dend, k = 3, order_clusters_as_data = FALSE) 

# Updated color mapping:
# cutree 1 = Cluster 2 (dark red)
# cutree 2 = Late Cluster 1 (light blue)
# cutree 3 = Early Cluster 1 (dark blue)
nm2 <- c("1" = "#B2182B", "2" = "#92C5DE", "3" = "#2166AC")

# export cluster year data for use in coloring other plots:
#saveRDS(k_3, "plots/Figure4_Subplots/k_3.RDS")

cols = nm2[as.character(k_3[Mesh_genus_yearly$year])]
cols

options(repr.plot.width=5, repr.plot.height=10)

DENDplot <- dend %>%
  set("labels_cex", 1.) %>%
  hang.dendrogram(hang_height = .1) %>% 
  set("leaves_pch", 19) %>% 
  set("leaves_col", cols, order_value = TRUE) %>% 
  set("leaves_bg", cols, order_value = TRUE) %>%
  set("branches_lwd", 0.5) 

DENDplot %>% plot(horiz=T)

DendroPLOT <- ggplot(DENDplot, horiz=FALSE, offset_labels=-0.03) + 
  theme_cowplot(font_size=20) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.clip = "off"
  ) + 
  scale_y_continuous(limits=c(0., 0.7)) + 
  scale_x_reverse()

DendroPLOT


# calculate Turnover
ds_genus_yearly$year <- as.integer(ds_genus_yearly$year)

# Calculate relative total turnover within replicates
total.res <- turnover(df=ds_genus_yearly,  
                      time.var = "year", 
                      species.var = "Genus",
                      abundance.var = "Total")

total.res$group = k_3[as.character(total.res$year)]
total.res$col = nm2[as.character(total.res$group)]
total.res$year

total.res$date <- as.Date(as.character(total.res$year), format="%Y")-135

total.res$col

datescale <- scale_x_date(date_labels = "%Y", expand=c(0,0), date_minor_breaks = "1 year", 
                          breaks=c(as.Date("1996/1/1"),as.Date("2000/1/1"),as.Date("2005/1/1"),as.Date("2010/1/1"),
                                   as.Date("2015/1/1"),as.Date("2017/1/1")), 
                          guide = guide_axis(minor.ticks = TRUE))

TurnoverPLOT <- ggplot(data=total.res, aes(x=date, y=total)) + 
  geom_line() + 
  geom_point(aes(col=col), size=3) + 
  theme_cowplot(font_size=20) +
  scale_colour_identity() +
  ylab("Species Turnover") + 
  xlab("Date [years]") +
  labs(colour="Group") + 
  guides(colour="none") + 
  datescale

TurnoverPLOT


### Extract dominant genus per Cluster from yearly aggregate community data

# Early Cluster 1: 1996-2003 (dark blue, cutree = 3)
GroupEarlyCluster1 <- ds_genus_yearly %>%  
  filter(year >= 1996 & year <= 2003)

# Late Cluster 1: 2014-2016 (light blue, cutree = 2)
GroupLateCluster1 <- ds_genus_yearly %>%  
  filter(year >= 2014 & year <= 2016)

# Cluster 2: 2004-2013 (dark red, cutree = 1)
GroupCluster2 <- ds_genus_yearly %>%  
  filter(year >= 2004 & year <= 2013)

returnTop5 <- function(x) {
  top5 <- x %>% group_by(Genus) %>% summarize(Full_sum = sum(Total)) %>% slice_max(order_by=Full_sum, n=5) 
  return(top5)
}

# Top 5 genera for each cluster
returnTop5(GroupEarlyCluster1)
returnTop5(GroupLateCluster1)
returnTop5(GroupCluster2)



################## Sup Plot - NMDS of Jaccard yearly ####################

# Run NMDS with Jaccard distance on binary (presence/absence) data
set.seed(42)
nmds_jac <- metaMDS(Mesh_genus_yearly[-1], distance = "jaccard", binary = TRUE, k = 2, trymax = 100)

# Check stress
print(paste("Stress:", round(nmds_jac$stress, 3)))

# Create scores dataframe
nmds_scores_jac <- as.data.frame(scores(nmds_jac, display = "sites"))
nmds_scores_jac$year <- as.integer(Mesh_genus_yearly$year)

# Map cluster assignments with descriptive labels
cluster_labels <- c("1" = "Cluster 2", "2" = "Late Cluster 1", "3" = "Early Cluster 1")
nmds_scores_jac$cluster <- factor(cluster_labels[as.character(k_3[as.character(nmds_scores_jac$year)])],
                                  levels = c("Early Cluster 1", "Cluster 2", "Late Cluster 1"))

# Define colors matching cluster labels
cluster_colors <- c("Early Cluster 1" = "#2166AC", 
                    "Cluster 2" = "#B2182B", 
                    "Late Cluster 1" = "#92C5DE")

# Create plot
NMDS_plot_jac <- ggplot(nmds_scores_jac, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = cluster), shape = 21, size = 4, color = "black", stroke = 0.5) +
  geom_text(aes(label = year), vjust = -1, size = 3.5) +
  scale_fill_manual(values = cluster_colors) +
  theme_cowplot(font_size = 14) +
  labs(title = paste("NMDS (Jaccard, presence/absence) — Stress:", round(nmds_jac$stress, 3)),
       fill = "Period") +
  theme(legend.position = "right")

NMDS_plot_jac