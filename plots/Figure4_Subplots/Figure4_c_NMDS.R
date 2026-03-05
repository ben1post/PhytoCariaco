library(tidyverse)
require(cowplot)
library(cowplot)
library(ggpubr)
library(vegan)
library(ggrepel)

# read interpolated phytoplankton cell count data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

unique(phyto_counts$FuncGroup)

# =============================================================================
# OPTION: Set to TRUE to exclude years without cluster assignment (1995, 2017)
# =============================================================================
exclude_unassigned <- TRUE

# =============================================================================
# PREPARE DATA FOR NMDS
# =============================================================================

# get Genus Matrix:
ds_genus <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, date) %>%
  summarise(Total = sum(counts), .groups = "drop") %>%
  arrange(desc(date))
tail(ds_genus, n=10)

ds_genus$year <- format(ds_genus$date, format="%Y")

ds_genus = ds_genus %>% filter(year>=1995 & year<=2017) %>% select(-year)

Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)

Mesh_genus_noemtpyrows = Mesh_genus[rowSums(Mesh_genus[, -1])>0, ]

# =============================================================================
# GET ENVIRONMENTAL DATA
# =============================================================================

CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")

CARIACO$u10_negative = -CARIACO$u10

new_Mesh_genus = Mesh_genus_noemtpyrows
new_Mesh_genus$time_month = format(new_Mesh_genus$date, format="%m-%Y")

CARIACO_dat_joined <- list(CARIACO, 
                           new_Mesh_genus
) %>% 
  reduce(full_join, by = c("time_month"))

sel_env_factors = c("u10_negative", "Salinity_bottles", "sst_10m", "Isotherm_21", "NO3_merged", "PO4_merged", "SiO4_merged", "AMO", "MEIv2", "tp")

firstSpec = as.numeric(which(names(CARIACO_dat_joined)=="Acanthoica"))
lastSpec = as.numeric(which(names(CARIACO_dat_joined)=="Zygosphaera"))
print(c(firstSpec, lastSpec))

Full_Matrix <- CARIACO_dat_joined %>%
  select("date", all_of(sel_env_factors), firstSpec:lastSpec)

Comp_Full_Matrix <- Full_Matrix[complete.cases(Full_Matrix),]

nfac = length(sel_env_factors) + 2
nfac

# =============================================================================
# VALIDATE GENUS COLUMN SELECTION
# =============================================================================

expected_last_genus_col <- ncol(Comp_Full_Matrix)
n_genus_cols <- expected_last_genus_col - nfac + 1

message(paste("Validated:", n_genus_cols, "genus columns selected (columns", nfac, "to", expected_last_genus_col, ")"))

# =============================================================================
# PREPARE ENVIRONMENTAL AND GENUS MATRICES
# =============================================================================

Env_Matrix <- Comp_Full_Matrix[, sel_env_factors]

names(Env_Matrix)[1] <- "Wind speed"
names(Env_Matrix)[2] <- "Salinity"
names(Env_Matrix)[3] <- "SST"
names(Env_Matrix)[4] <- "Isotherm Depth"
names(Env_Matrix)[5] <- "NO3"
names(Env_Matrix)[6] <- "PO4"
names(Env_Matrix)[7] <- "SiO4"
names(Env_Matrix)[8] <- "AMO"
names(Env_Matrix)[9] <- "MEIv2"
names(Env_Matrix)[10] <- "Precipitation"

# Select all genus columns dynamically (no hardcoded offset)
Genus_Matrix <- Comp_Full_Matrix[, nfac:ncol(Comp_Full_Matrix)]

# Remove rare genera (present in <= 5 samples)
Genus_Matrix_2 <- Genus_Matrix %>% select_if(colSums(. > 0) > 5)
message(paste("Genera before filtering:", ncol(Genus_Matrix), "| After filtering:", ncol(Genus_Matrix_2)))

# =============================================================================
# NMDS ANALYSIS (Bray-Curtis, cube-root transformed)
# =============================================================================

# Cube-root transformation to down-weight dominant taxa
m_com = as.matrix(Genus_Matrix_2^(1/3))

# NMDS with Bray-Curtis distance, k=2 dimensions
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", k=2, trymax=200)
nmds

# =============================================================================
# ENVIRONMENTAL VECTOR FITTING
# =============================================================================

en = envfit(nmds, Env_Matrix, permutations = 999, na.rm = TRUE)
en

plot(nmds)
plot(en)

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

# Extract NMDS scores (x and y coordinates) for sites
data.scores_2 = as.data.frame(scores(nmds)$sites)

# Add year column
data.scores_2$Year = format(Comp_Full_Matrix$date, "%Y")

# Load cluster assignments
k_3 <- readRDS("plots/Figure4_Subplots/k_3.RDS")

# Assign clusters, with "unassigned" for years outside 1996-2016
data.scores_2$cluster <- sapply(data.scores_2$Year, function(y) {
  if (y %in% names(k_3)) {
    as.character(k_3[y])
  } else {
    "unassigned"
  }
})

head(data.scores_2)

# Optionally exclude unassigned years
if (exclude_unassigned) {
  data.scores_2 <- data.scores_2 %>% filter(cluster != "unassigned")
  message("Excluded unassigned years (1995, 2017) from plot")
} else {
  message("Including unassigned years (1995, 2017) in grey")
}

# =============================================================================
# NMDS PLOT
# =============================================================================

options(repr.plot.width=10, repr.plot.height=10)

# Color vector for 3 clusters + unassigned
# Cluster 2 (2004-2013): dark red
# Late Cluster 1 (2014-2016): light blue
# Early Cluster 1 (1996-2003): dark blue
# Unassigned (1995, 2017): grey
colsvec <- c("1" = "#B2182B", "2" = "#92C5DE", "3" = "#2166AC", "unassigned" = "grey50")


NMDSplot <- ggplot(data.scores_2, aes(x = NMDS1, y = NMDS2)) + 
  
  geom_point(size = 2, aes(colour = cluster)) + 
  
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  
  labs(x = "NMDS1", colour = "Cluster", y = "NMDS2", shape = "Group") + 
  
  scale_colour_manual(values = colsvec) + 
  
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2, label = row.names(en_coord_cont)), 
            colour = "grey30", fontface = "bold", size = 4.5, hjust = -0.1, vjust = 0.5) +
  
  theme_cowplot(font_size = 20) + 
  
  guides(colour = "none")

NMDSplot
