library(tidyverse)

# Read Phytoplankton Data:
phyto_counts_old <- readRDS("data/processed/PhytoplanktonInterpolatedCounts_preReview.RDS")
phyto_counts_new <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

str(phyto_counts_old)
str(phyto_counts_new)

# Read Combined Data
CARIACO_old <- readRDS("data/processed/CARIACO_EnvData_combined_preReview.rds")
CARIACO_new <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


str(CARIACO_old)
str(CARIACO_new)




# =============================================================================
# DIAGNOSE JUNE 2016 NAs
# =============================================================================

# Read raw data
ds <- read.csv("data/BCO-DMO/phytoplankton.csv")
ds$date <- lubridate::parse_date_time(ds$Datetime_UTC, orders="%Y-%m-%d[.]H:M", tz='UTC')

# =============================================================================
# 1. Is June 2016 in the raw data at all?
# =============================================================================

cat("=== Q1: Is June 2016 in raw data? ===\n")
june_2016_raw <- ds %>% filter(format(date, "%Y-%m") == "2016-06")
cat("Rows in June 2016:", nrow(june_2016_raw), "\n")
cat("Unique dates:", as.character(unique(june_2016_raw$date)), "\n")
cat("Unique AphiaIDs:", n_distinct(june_2016_raw$AphiaID), "\n")

# =============================================================================
# 2. What do the depth values look like for June 2016?
# =============================================================================

cat("\n=== Q2: Depth data quality for June 2016 ===\n")
depth_cols <- c("d_1m", "d_7m", "d_15m", "d_25m", "d_35m", "d_55m", "d_75m", "d_100m")

# Convert to numeric
for (col in depth_cols) {
  june_2016_raw[[col]] <- as.numeric(june_2016_raw[[col]])
}

# Count NAs vs 0s vs positive values per depth
cat("Per depth column:\n")
for (col in depth_cols) {
  n_na <- sum(is.na(june_2016_raw[[col]]))
  n_zero <- sum(june_2016_raw[[col]] == 0, na.rm = TRUE)
  n_pos <- sum(june_2016_raw[[col]] > 0, na.rm = TRUE)
  cat(sprintf("  %s: %d NA, %d zero, %d positive\n", col, n_na, n_zero, n_pos))
}

# =============================================================================
# 3. Compare to a "normal" cruise
# =============================================================================

cat("\n=== Q3: Compare to Jan 2016 (normal cruise) ===\n")
jan_2016_raw <- ds %>% filter(format(date, "%Y-%m") == "2016-01")
cat("Jan 2016 rows:", nrow(jan_2016_raw), "\n")

for (col in depth_cols) {
  jan_2016_raw[[col]] <- as.numeric(jan_2016_raw[[col]])
}

cat("Per depth column (Jan 2016):\n")
for (col in depth_cols) {
  n_na <- sum(is.na(jan_2016_raw[[col]]))
  n_zero <- sum(jan_2016_raw[[col]] == 0, na.rm = TRUE)
  n_pos <- sum(jan_2016_raw[[col]] > 0, na.rm = TRUE)
  cat(sprintf("  %s: %d NA, %d zero, %d positive\n", col, n_na, n_zero, n_pos))
}

# =============================================================================
# 4. Check if June 2016 taxa are in the corrected names table
# =============================================================================

cat("\n=== Q4: Are June 2016 taxa in corrected names table? ===\n")
occurrence_corrected <- read.csv("data/processed/corrected_phyto_names_and_ids.csv")

june_aphia <- unique(june_2016_raw$AphiaID)
corrected_aphia <- unique(occurrence_corrected$AphiaID)

cat("June 2016 unique AphiaIDs:", length(june_aphia), "\n")
cat("AphiaIDs in corrected table:", length(corrected_aphia), "\n")
cat("June 2016 AphiaIDs NOT in corrected table:", 
    length(setdiff(june_aphia, corrected_aphia)), "\n")

missing_aphia <- setdiff(june_aphia, corrected_aphia)
if (length(missing_aphia) > 0 & length(missing_aphia) <= 20) {
  cat("Missing AphiaIDs:\n")
  print(missing_aphia)
}

# =============================================================================
# 5. Trace through the pipeline for one specific NA taxon
# =============================================================================

cat("\n=== Q5: Trace a specific NA case ===\n")

# Pick a taxon that has NA in June 2016
na_taxon <- phyto_counts_new %>%
  filter(format(date, "%Y-%m") == "2016-06" & is.na(counts)) %>%
  slice(1)

cat("Example NA taxon:\n")
print(na_taxon)

example_aphia <- na_taxon$AphiaID

# Is this taxon in June 2016 raw data?
cat("\nIs AphiaID", example_aphia, "in June 2016 raw data?\n")
example_in_june <- june_2016_raw %>% filter(AphiaID == example_aphia)
cat("Rows:", nrow(example_in_june), "\n")

if (nrow(example_in_june) > 0) {
  cat("Depth values:\n")
  print(example_in_june[, depth_cols])
}

# Is this taxon in the corrected names?
cat("\nIs AphiaID", example_aphia, "in corrected names table?\n")
example_corrected <- occurrence_corrected %>% filter(AphiaID == example_aphia)
cat("Rows:", nrow(example_corrected), "\n")
if (nrow(example_corrected) > 0) {
  print(example_corrected[, c("AphiaID", "ScientificName_corrected", "CorrectedAphiaID")])
}

# =============================================================================
# 6. Check what happens after the merge in the processing script
# =============================================================================

cat("\n=== Q6: Check merge results ===\n")

# Recreate the key merge step from your script
AphiaIDcorrected <- data.frame(
  "AphiaID" = as.integer(occurrence_corrected$AphiaID), 
  "CorrectedAphiaID" = as.integer(occurrence_corrected$CorrectedAphiaID), 
  "SpeciesNameCleaned" = occurrence_corrected$SpeciesNameCleaned
)

AphiaIDcorrected_rmdp <- AphiaIDcorrected %>% 
  distinct(AphiaID, SpeciesNameCleaned, .keep_all = TRUE)

# How many June 2016 taxa survive the merge?
june_merged <- merge(june_2016_raw, AphiaIDcorrected_rmdp, 
                     by = c("AphiaID", "SpeciesNameCleaned"), all = FALSE)
cat("June 2016 rows after merge with corrected names:", nrow(june_merged), "\n")
cat("June 2016 rows LOST in merge:", nrow(june_2016_raw) - nrow(june_merged), "\n")

# =============================================================================
# 7. Final check: Count of taxa with values vs NAs in new data for June 2016
# =============================================================================

cat("\n=== Q7: Summary of June 2016 in final output ===\n")
june_final <- phyto_counts_new %>% filter(format(date, "%Y-%m") == "2016-06")
cat("Total rows:", nrow(june_final), "\n")
cat("Rows with counts (not NA):", sum(!is.na(june_final$counts)), "\n")
cat("Rows with NA:", sum(is.na(june_final$counts)), "\n")
cat("Rows with counts = 0:", sum(june_final$counts == 0, na.rm = TRUE), "\n")
cat("Rows with counts > 0:", sum(june_final$counts > 0, na.rm = TRUE), "\n")
  
  
  
  
  