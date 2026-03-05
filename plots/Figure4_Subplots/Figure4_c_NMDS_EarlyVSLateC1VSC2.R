library(tidyverse)
library(vegan)
library(cowplot)
library(ggrepel)

# --- 1. LOAD AND PREPARE DATA (same as your scripts) ---
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

ds_genus <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, date) %>%
  summarise(Total = sum(counts), .groups = "drop") %>%
  arrange(desc(date))

ds_genus$year <- format(ds_genus$date, format = "%Y")
ds_genus <- ds_genus %>% select(-year)#ds_genus %>% filter(year >= 1996 & year <= 2016) %>% select(-year)

Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, 
                          values_fn = sum, values_fill = 0)
Mesh_genus_noemtpyrows <- Mesh_genus[rowSums(Mesh_genus[, -1]) > 0, ]

# --- 2. JOIN WITH ENVIRONMENTAL DATA ---
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")
CARIACO$u10_negative <- -CARIACO$u10

new_Mesh_genus <- Mesh_genus_noemtpyrows
new_Mesh_genus$time_month <- format(new_Mesh_genus$date, format = "%m-%Y")

CARIACO_dat_joined <- list(CARIACO, new_Mesh_genus) %>% 
  reduce(full_join, by = "time_month")

sel_env_factors <- c("u10_negative", "Salinity_bottles", "sst_10m", "Isotherm_21", 
                     "NO3_merged", "PO4_merged", "SiO4_merged", "AMO", "MEIv2", "tp")

firstSpec <- as.numeric(which(names(CARIACO_dat_joined) == "Acanthoica"))
lastSpec <- as.numeric(which(names(CARIACO_dat_joined) == "Zygosphaera"))

Full_Matrix <- CARIACO_dat_joined %>%
  select("date", all_of(sel_env_factors), firstSpec:lastSpec)

Comp_Full_Matrix <- Full_Matrix[complete.cases(Full_Matrix), ]

nfac <- length(sel_env_factors) + 2

# --- 3. PREPARE MATRICES ---
Env_Matrix <- Comp_Full_Matrix[, sel_env_factors]
names(Env_Matrix) <- c("Wind speed", "Salinity", "SST", "Isotherm Depth", 
                       "NO3", "PO4", "SiO4", "AMO", "MEIv2", "Precipitation")

Genus_Matrix <- Comp_Full_Matrix[, nfac:(nfac + 222)]
Genus_Matrix_2 <- Genus_Matrix %>% select_if(colSums(. > 0) > 5)

# Cube root transformation
m_com <- as.matrix(Genus_Matrix_2^(1/3))

# --- 4. RUN NMDS ---
set.seed(123)
nmds <- metaMDS(m_com, distance = "bray", k = 2, trymax = 200)

#nmds = metaMDS(m_com, distance = "jaccard", binary = TRUE, k=3, trymax=200)

# --- 5. ENVFIT ---
en <- envfit(nmds, Env_Matrix, permutations = 999, na.rm = TRUE)
en_coord_cont <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

# --- 6. EXTRACT SCORES AND ASSIGN PERIODS ---
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$Year <- format(Comp_Full_Matrix$date, "%Y")
data.scores$date <- Comp_Full_Matrix$date

# Load cluster assignments
k_3 <- readRDS("plots/Figure4_Subplots/k_3.RDS")
data.scores$cluster <- as.character(k_3[data.scores$Year])

# Create three-period classification
data.scores <- data.scores %>%
  mutate(
    Period = case_when(
      cluster == "3" & as.numeric(Year) <= 2003 ~ "Early Cluster 1 (1996-2003)",
      cluster == "1" ~ "Cluster 2 (2004-2013)",
      cluster == "2" & as.numeric(Year) >= 2014 ~ "Late Cluster 1 (2014-2016)",
      TRUE ~ "Other"
    ),
    Period = factor(Period, levels = c("Early Cluster 1 (1996-2003)", 
                                       "Cluster 2 (2004-2013)", 
                                       "Late Cluster 1 (2014-2016)"))
  )

# Check assignment
cat("=== Period assignments ===\n")
print(table(data.scores$Period))
cat("\n=== Years in each period ===\n")
data.scores %>% group_by(Period) %>% summarise(Years = paste(sort(unique(Year)), collapse = ", ")) %>% print()

# --- 7. DEFINE COLORS ---
# Blue shades for Cluster 1 (early = darker, late = lighter), Red for Cluster 2
period_cols <- c(
  "Early Cluster 1 (1996-2003)" = "#2166AC",   # Dark blue
  "Cluster 2 (2004-2013)" = "#B2182B",          # Red
  "Late Cluster 1 (2014-2016)" = "#92C5DE"      # Light blue
)

# --- 8. CREATE PLOT ---
NMDSplot_periods <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  
  # Points colored by period
  geom_point(aes(colour = Period), size = 2.5, alpha = 0.8) + 
  scale_colour_manual(values = period_cols) +
  
  # Environmental vectors
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
                  colour = "grey30", fontface = "bold", 
                  label = row.names(en_coord_cont), size = 4.5) +
  
  # Labels and theme
  labs(x = "NMDS1", y = "NMDS2", colour = "Period") +
  theme_cowplot(font_size = 20) +
  theme(
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.position = "right"
  )

print(NMDSplot_periods)

# --- 9. OPTIONAL: Add convex hulls around periods ---
# This visually emphasizes the groupings

# Calculate hulls
hulls <- data.scores %>%
  group_by(Period) %>%
  slice(chull(NMDS1, NMDS2))

NMDSplot_periods_hulls <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  
  # Convex hulls (semi-transparent)
  geom_polygon(data = hulls, aes(fill = Period), alpha = 0.2) +
  scale_fill_manual(values = period_cols) +
  
  # Points
  geom_point(aes(colour = Period), size = 2.5) + 
  scale_colour_manual(values = period_cols) +
  
  # Environmental vectors
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
                  colour = "grey30", fontface = "bold", 
                  label = row.names(en_coord_cont), size = 4.5) +
  
  labs(x = "NMDS1", y = "NMDS2", colour = "Period", fill = "Period") +
  theme_cowplot(font_size = 20) +
  theme(
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.position = "right"
  )

print(NMDSplot_periods_hulls)

# --- 10. QUANTIFY OVERLAP ---
# Calculate centroid distances between periods
cat("\n=== Period centroids in NMDS space ===\n")
centroids <- data.scores %>%
  group_by(Period) %>%
  summarise(
    NMDS1_mean = mean(NMDS1),
    NMDS2_mean = mean(NMDS2),
    n = n(),
    .groups = "drop"
  )
print(centroids)

# Distance between centroids
cat("\n=== Euclidean distances between period centroids ===\n")
cat(sprintf("Early C1 to Cluster 2: %.3f\n", 
            sqrt((centroids$NMDS1_mean[1] - centroids$NMDS2_mean[2])^2 + 
                   (centroids$NMDS2_mean[1] - centroids$NMDS2_mean[2])^2)))
cat(sprintf("Late C1 to Cluster 2: %.3f\n", 
            sqrt((centroids$NMDS1_mean[3] - centroids$NMDS1_mean[2])^2 + 
                   (centroids$NMDS2_mean[3] - centroids$NMDS2_mean[2])^2)))
cat(sprintf("Early C1 to Late C1: %.3f\n", 
            sqrt((centroids$NMDS1_mean[1] - centroids$NMDS1_mean[3])^2 + 
                   (centroids$NMDS2_mean[1] - centroids$NMDS2_mean[3])^2)))






#########################

# --- 1. SUBSET DATA ---
# We only want to compare Early vs Late Cluster 1
target_periods <- c("Early Cluster 1 (1996-2003)", "Late Cluster 1 (2014-2016)")

# Filter the community matrix (raw transformed counts) and the grouping factor
subset_indices <- data.scores$Period %in% target_periods
m_com_subset <- m_com[subset_indices, ]
groups_subset <- factor(data.scores$Period[subset_indices])

# --- 2. PERMANOVA (Location Test) ---
# Testing H0: Centroids of Early and Late periods are identical
set.seed(123)
permanova_res <- adonis2(m_com_subset ~ groups_subset, 
                         method = "bray", 
                         permutations = 999)

cat("=== PERMANOVA Results (Early vs Late) ===\n")
print(permanova_res)

# --- 3. HOMOGENEITY OF DISPERSION (Assumption Check) ---
# Testing H0: Variances (spread) of Early and Late periods are identical
dist_matrix <- vegdist(m_com_subset, method = "bray")
dispersion <- betadisper(dist_matrix, groups_subset)
dispersion_test <- permutest(dispersion, pairwise = TRUE, permutations = 999)

cat("\n=== Dispersion Test Results ===\n")
print(dispersion_test)

# Optional: Visualize dispersion
# plot(dispersion, main = "Multivariate Dispersion: Early vs Late")


###### --- MODIFIED PLOT: Late Cluster 1 points colored by year ---######

# Create a combined color/shape scheme
data.scores <- data.scores %>%
  mutate(
    # For coloring: use Period for Early/C2, use Year for Late C1
    PlotColor = case_when(
      Period == "Early Cluster 1 (1996-2003)" ~ "Early C1",
      Period == "Cluster 2 (2004-2013)" ~ "Cluster 2",
      Period == "Late Cluster 1 (2014-2016)" ~ Year,  # Individual years
      TRUE ~ "Other"
    ),
    PlotColor = factor(PlotColor, levels = c("Early C1", "Cluster 2", 
                                             "2014", "2015", "2016"))
  )

# High-contrast color palette
detailed_cols <- c(
  "Early C1" = "#2166AC",      # Dark blue
  "Cluster 2" = "#B2182B",     # Red
  "2014" = "#FF7F00",          # Orange
  "2015" = "#6A3D9A",          # Purple
  "2016" = "#33A02C"           # Green
)

# Create the plot
NMDSplot_late_years <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  
  # Points
  geom_point(aes(colour = PlotColor), size = 2.5, alpha = 0.8) + 
  scale_colour_manual(
    values = detailed_cols,
    labels = c("Early C1 (1996-2003)", "Cluster 2 (2004-2013)", 
               "2014", "2015", "2016")
  ) +
  
  # Environmental vectors
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
                  colour = "grey30", fontface = "bold", 
                  label = row.names(en_coord_cont), size = 4.5,
                  max.overlaps = 15) +
  
  labs(x = "NMDS1", y = "NMDS2", colour = "Period / Year") +
  theme_cowplot(font_size = 20) +
  theme(
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.position = "right"
  )

print(NMDSplot_late_years)

# --- Calculate year-specific centroids for 2014-2016 ---
cat("\n=== Late period centroids by year ===\n")
late_centroids <- data.scores %>%
  filter(Period == "Late Cluster 1 (2014-2016)") %>%
  group_by(Year) %>%
  summarise(
    NMDS1_mean = mean(NMDS1),
    NMDS2_mean = mean(NMDS2),
    n = n(),
    .groups = "drop"
  )
print(late_centroids)

# Add overall centroids for comparison
cat("\n=== Distance of each late year to Early C1 vs Cluster 2 centroids ===\n")
early_centroid <- c(centroids$NMDS1_mean[1], centroids$NMDS2_mean[1])
c2_centroid <- c(centroids$NMDS1_mean[2], centroids$NMDS2_mean[2])

late_centroids <- late_centroids %>%
  mutate(
    dist_to_early = sqrt((NMDS1_mean - early_centroid[1])^2 + 
                           (NMDS2_mean - early_centroid[2])^2),
    dist_to_c2 = sqrt((NMDS1_mean - c2_centroid[1])^2 + 
                        (NMDS2_mean - c2_centroid[2])^2),
    closer_to = ifelse(dist_to_early < dist_to_c2, "Early C1", "Cluster 2")
  )
print(late_centroids)

