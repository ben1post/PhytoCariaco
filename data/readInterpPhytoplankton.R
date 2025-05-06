library(tidyverse)

# import interpolation functions
source('data/interpolateData.r')

# Read phytoplankton count data
ds <- read.csv("data/BCO-DMO/phytoplankton.csv")  #, na.strings=c("nd","?"))

# Fix date, Aphia ID and depth columns
ds$date = parse_date_time(ds$Datetime_UTC, orders="%Y-%m-%d[.]H:M", tz='UTC')
ds$d_7m = as.numeric(ds$d_7m)
ds$d_75m = as.numeric(ds$d_75m)
ds$d_100m = as.numeric(ds$d_100m)
ds$AphiaID = as.integer(ds$AphiaID)

# read annotated and corrected pyhytoplankton data (table with manually corrected names and IDs)
occurrence_corrected <- read.csv("data/processed/corrected_phyto_names_and_ids.csv", stringsAsFactor=FALSE)

# trim dataframe to relevant columns
AphiaIDcorrected <- data.frame("AphiaID" = as.integer(occurrence_corrected$AphiaID), 
                               "TaxonRank"=as.character(occurrence_corrected$IdentifiedRank),
                               "CorrectedAphiaID"=as.integer(occurrence_corrected$CorrectedAphiaID), 
                               "SpeciesNameCleaned"=occurrence_corrected$SpeciesNameCleaned, 
                               "ScientificName_corrected"=occurrence_corrected$ScientificName_corrected,
                               "FuncGroup"=occurrence_corrected$FuncGroup,
                               "Genus"=occurrence_corrected$Genus, stringsAsFactors=FALSE)


# Filter out duplicates
AphiaIDcorrected_rmdp <- AphiaIDcorrected %>% distinct(AphiaID, ScientificName_corrected, .keep_all = TRUE)

# Merge the original data and the corrected by species name and Aphia ID
ds_phytoMergedCorrected <- merge(ds, AphiaIDcorrected_rmdp, by=c("AphiaID", "SpeciesNameCleaned"), all=TRUE)

# add nanoflagellate data with placeholder AphiaID:
ds_phytoMergedCorrected[ds_phytoMergedCorrected$SpeciesNameCleaned=="nanoflagellates",]$AphiaID <- 1

# replace Aphia IDs from original dataset with corrected annotated IDs, then remove NA Aphia IDs (only identified samples)
ds_FG <- ds_phytoMergedCorrected %>% mutate(AphiaID = coalesce(AphiaID, CorrectedAphiaID)) %>% drop_na(AphiaID)

# extract metadata of unique identified units for later merging with interpolated counts via Aphia ID:
ds_FG_uniqueInfo <- ds_FG  %>% group_by(AphiaID) %>% 
  summarize(TaxonRank=first(TaxonRank, na_rm = TRUE), FuncGroup=first(FuncGroup, na_rm = TRUE), 
            Genus=first(Genus, na_rm = TRUE), ScientificName_corrected=first(ScientificName_corrected, na_rm = TRUE))


# prepare data for interpolation:
ds_phyInt <- rbind(data.frame(val=ds_FG$d_1m, depth=1, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_7m, depth=7, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_15m, depth=15, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_25m, depth=25, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_55m, depth=55, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_75m, depth=75, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank),
                   data.frame(val=ds_FG$d_100m, depth=100, date=ds_FG$date, AphiaID=ds_FG$AphiaID, Genus=ds_FG$Genus, TaxonRank=ds_FG$TaxonRank))


# Pivot dataframe to have AphiaID as columns - Missing Data is set to 0 (none observed)
Mesh_phyInt <- pivot_wider(ds_phyInt, names_from = AphiaID, values_from = val, values_fn = sum, values_fill = 0)

# Pivot dataframe to have Genera as columns - Missing Data is set to 0 (none observed)
ds_genus <- ds_phyInt %>% filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% group_by(Genus, date, depth) %>% summarise(Total = sum(val)) %>% arrange(date)
Mesh_genus_phyInt <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)

# Export raw counts per Genus for depth-interval interpolation in "calculatePhytoDepthIntervals.R":
saveRDS(Mesh_genus_phyInt, "data/processed/PhytoplanktonRawGenusCounts.RDS")

# Get unique Aphia IDs to iterate over
phyto_AphiaID = na.omit(as.character(unique(ds_phyInt$AphiaID)))

# Interpolate
getPhytoInterpCounts <- function(depth_from=0, depth_to=100, noofNA=20){
  phyto_temp_store = list()
  
  for (variable in phyto_AphiaID) {
    # interpolation algorithm: oce-rr
    phyto_temp_store[[variable]] <- interpolateData(Mesh_phyInt, variable, depth_from=depth_from, depth_to=depth_to, noofNA=noofNA, int_func='unesco')
    names(phyto_temp_store[[variable]])[1] <- variable
  }
  
  phyto_ds_cleaned <- phyto_temp_store %>% 
    reduce(full_join, by = "date") %>% na.omit()
  
  phyto_ds_cleaned_pivot <- phyto_ds_cleaned %>% 
    pivot_longer(cols=-date, names_to = "AphiaID", values_to = "counts") 
  
  phyto_ds_cleaned_pivot$AphiaID = as.integer(phyto_ds_cleaned_pivot$AphiaID)
  
  ds_mergedIntCounts <- merge(phyto_ds_cleaned_pivot, ds_FG_uniqueInfo, by=c("AphiaID"), all=TRUE) %>% 
    select(date, AphiaID, counts, ScientificName_corrected, Genus, FuncGroup, TaxonRank)
  
  return(ds_mergedIntCounts)    
}

ds_mergedIntCounts <- getPhytoInterpCounts()

str(ds_mergedIntCounts)

# Save interpolated phytoplankton counts
saveRDS(ds_mergedIntCounts, "data/processed/PhytoplanktonInterpolatedCounts.RDS")
