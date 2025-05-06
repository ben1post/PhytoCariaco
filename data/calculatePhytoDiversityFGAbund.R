library(tidyverse)
library(vegan)

# Read interpolated counts (from script "readInterpPhytoplankton.R")
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
# Remove NAs
phyto_counts = phyto_counts[complete.cases(phyto_counts$counts),]
str(phyto_counts)
# Fix date format
phyto_counts$realdate <- phyto_counts$date 
phyto_counts$date <- as.Date(format(phyto_counts$date,"%Y-%m-15"))

unique(phyto_counts$ScientificName_corrected[phyto_counts$TaxonRank == "Group"])
unique(phyto_counts$FuncGroup)

# get Genus Matrix:
ds_genus <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, date) %>%
  summarise(Total = sum(counts))  %>%
  arrange(desc(date))
#tail(ds_genus, n=10)

# Pivot Genus Matrix for calculating diversity indices (Aphia ID as columns)
Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)

# Calculate Diversity indices
Shannon_gen <- diversity(Mesh_genus[-1])
Pielou_gen <- Shannon_gen/log(specnumber(Mesh_genus[-1]))
GenusRichness <- apply(Mesh_genus[-1]>0,1,sum)


# Get Genus Matrix incl. all group counts (also nanoflagellate counts) 
ds_genus_full <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species" | TaxonRank == "Group") %>%
  group_by(Genus, date, FuncGroup) %>%
  summarise(Total = sum(counts))  %>%
  arrange(desc(date))

# get abundance for each date for each functional group
Gen_abund_FG <- ds_genus_full %>% group_by(date, FuncGroup) %>% 
  summarise(abundance = sum(Total))  %>%
  arrange(desc(date))

# Pivot data frame for further processing
Gen_abund_FG_wide <- pivot_wider(Gen_abund_FG, names_from = FuncGroup, values_from = abundance)


# Combine data into single data frame
diversity_ds = data.frame("time_month"=format(Mesh_genus[[1]], format="%m-%Y"), GenusRichness, Shannon_gen, Pielou_gen,
                          "Abundance_Dino"=Gen_abund_FG_wide$Dinoflagellata,
                          "Abundance_Diatom"=Gen_abund_FG_wide$Bacillariophyceae,
                          "Abundance_Hapto"=Gen_abund_FG_wide$Haptophyta,
                          "Abundance_Cyano"=Gen_abund_FG_wide$Cyanobacteria,
                          "Abundance_Crypto"=Gen_abund_FG_wide$Cryptophyte,
                          "Abundance_Chloro"=Gen_abund_FG_wide$Chlorophyte,
                          "Abundance_Nanoflagellate"=Gen_abund_FG_wide$Nanoflagellates)

# Export and save
saveRDS(diversity_ds, "data/processed/PhytoplanktonAbundanceDiversity.RDS")
