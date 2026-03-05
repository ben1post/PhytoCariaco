library(tidyverse)
library(vegan)

# Read interpolated counts
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

# Remove NAs
phyto_counts = phyto_counts[complete.cases(phyto_counts$counts),]

# Fix date format
phyto_counts$realdate <- phyto_counts$date 
phyto_counts$date <- as.Date(format(phyto_counts$date,"%Y-%m-15"))



#check for most abundant species:
# --- 3. MOST ABUNDANT SPECIES ---
# Filter for species-level records and sum across entire time series
species_abundance <- phyto_counts %>%
  filter(TaxonRank == "Species") %>%
  group_by(ScientificName_corrected) %>%
  summarise(
    Total_Counts = sum(counts),
    Mean_Counts = mean(counts),
    N_Observations = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Counts))

# View top 20 most abundant species
print(species_abundance, n = 20)

# Quick visualization of top 15
species_abundance %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(ScientificName_corrected, Total_Counts), y = Total_Counts, fill = FuncGroup)) +
  geom_col() +
  coord_flip() +
  labs(x = "Species", 
       y = "Total Cell Counts",
       title = "Top 15 Most Abundant Phytoplankton Species",
       fill = "Functional Group") +
  theme_minimal()

# --- 1. CALCULATE DIVERSITY ---

# Get Genus Matrix:
ds_genus <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>%
  group_by(Genus, date) %>%
  summarise(Total = sum(counts), .groups="drop")  %>% # ungroup to be safe
  arrange(desc(date))

# Pivot Genus Matrix
Mesh_genus <- pivot_wider(ds_genus, 
                          names_from = Genus, 
                          values_from = Total, 
                          values_fn = sum, 
                          values_fill = 0)

# Calculate Indices
# Safer to explicitly exclude the date column by name rather than index [-1]
count_matrix <- Mesh_genus %>% select(-date)

Shannon_gen <- diversity(count_matrix)
Pielou_gen  <- Shannon_gen / log(specnumber(count_matrix))
GenusRichness <- apply(count_matrix > 0, 1, sum)

# Create a safe Diversity Dataframe
df_diversity <- data.frame(
  date = Mesh_genus$date,
  GenusRichness = GenusRichness,
  Shannon_gen = Shannon_gen,
  Pielou_gen = Pielou_gen
)

# --- 2. CALCULATE ABUNDANCES ---

# Get Genus Matrix incl. all group counts
ds_genus_full <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species" | TaxonRank == "Group"| TaxonRank == "Phylum") %>%
  group_by(Genus, date, FuncGroup) %>%
  summarise(Total = sum(counts), .groups="drop") %>%
  arrange(desc(date))

# Get abundance for each date for each functional group
Gen_abund_FG <- ds_genus_full %>% 
  group_by(date, FuncGroup) %>% 
  summarise(abundance = sum(Total), .groups="drop")

# Pivot data frame
Gen_abund_FG_wide <- pivot_wider(Gen_abund_FG, 
                                 names_from = FuncGroup, 
                                 values_from = abundance,
                                 values_fill = 0) # Fill 0 if a group is missing on a date

# --- 3. SAFE MERGE (THE FIX) ---

# Join by DATE to ensure alignment, even if row counts differ
diversity_ds <- full_join(df_diversity, Gen_abund_FG_wide, by = "date") %>%
  mutate(time_month = format(date, format="%m-%Y")) %>%
  
  # Rename columns to match your desired output
  rename(
    Abundance_Dino = Dinoflagellata,
    Abundance_Diatom = Bacillariophyceae,
    Abundance_Hapto = Haptophyta,
    Abundance_Cyano = Cyanobacteria,
    Abundance_Crypto = Cryptophyte,
    Abundance_Chloro = Chlorophyte,
    Abundance_Nanoflagellate = Nanoflagellates
  ) %>%
  
  # Reorder columns (optional, puts time_month first)
  select(time_month, GenusRichness, Shannon_gen, Pielou_gen, starts_with("Abundance_"))

# Check structure
str(diversity_ds)

# Export and save
saveRDS(diversity_ds, "data/processed/PhytoplanktonAbundanceDiversity.RDS")
