library(tidyr)
library(dplyr)
library(lubridate)


# read combined dataset:
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")
CARIACO$year <- format(CARIACO$date, format="%Y")


# reformat and scale data variables for visualisation
CARIACO$u10_negative = -CARIACO$u10
#CARIACO$Isotherm_21_negative = -CARIACO$Isotherm_21

CARIACO$Chlorophyll_log = log1p(CARIACO$Chlorophyll)
CARIACO$Abundance_Diatom_log = log1p(CARIACO$Abundance_Diatom)
CARIACO$Abundance_Hapto_log = log1p(CARIACO$Abundance_Hapto)
CARIACO$Abundance_Dino_log = log1p(CARIACO$Abundance_Dino)
CARIACO$Abundance_Cyano_log = log1p(CARIACO$Abundance_Cyano)
CARIACO$Abundance_Nanoflagellate_log = log1p(CARIACO$Abundance_Nanoflagellate)

#### Diagnose Missing Months from Time series ####


# 1. Define the full timeline (1995-2017) ensuring we match your day-of-month format
# Your data uses the 15th of each month, so we match that to ensure joining works.
full_dates <- seq(as.Date("1995-01-15"), as.Date("2017-12-15"), by = "month")

# 2. Map variable codes to the readable names you used in the plot
# (We create a little lookup table from your existing vectors)
var_lookup <- data.frame(
  var_code = vars, 
  var_name = renamed_vars, 
  stringsAsFactors = FALSE
)

# 3. Create the robust summary
missing_table_robust <- CARIACO %>%
  # Keep only the date and the variables of interest
  select(date, all_of(vars)) %>%
  
  # Reshape to "Long" format so we can process all variables at once
  pivot_longer(cols = -date, names_to = "var_code", values_to = "value") %>%
  
  # --- THE CRITICAL STEP ---
  # This forces the dataframe to have a row for every combination of 
  # variable and date. If a month was missing, it is created here with value = NA.
  complete(var_code, date = full_dates) %>% 
  
  # Extract Year for grouping
  mutate(Year = format(date, "%Y")) %>%
  
  # Now filter to your desired range (just in case)
  filter(Year >= 1995 & Year <= 2017) %>%
  
  # Count NAs (This now counts both original NAs AND missing months)
  group_by(var_code, Year) %>%
  summarize(n_na = sum(is.na(value)), .groups = "drop") %>%
  
  # Join with your readable names
  left_join(var_lookup, by = "var_code") %>%
  
  # Apply the sorting order from your plot
  mutate(var_name = factor(var_name, levels = var_reorder)) %>%
  arrange(var_name) %>%
  
  # Clean up columns for final table
  select(var_name, Year, n_na) %>%
  
  # Pivot back to Wide format for the readable table
  pivot_wider(names_from = Year, values_from = n_na)

# 4. View and Export
print(missing_table_robust, n=21)

# 3. Helper function to format consecutive numbers as ranges
format_ranges <- function(x) {
  if (length(x) == 0) return("-")
  x <- sort(unique(x)) # Ensure sorted
  
  # Identify breaks in sequence to form groups
  # diff(x) == 1 means consecutive. cumulative sum of (diff != 1) gives group IDs
  groups <- cumsum(c(1, diff(x) != 1))
  
  # Split x into list of consecutive vectors
  split_x <- split(x, groups)
  
  # Format each group
  ranges <- sapply(split_x, function(v) {
    if (length(v) > 2) {
      # If 3 or more consecutive months, use Range format (Start-End)
      paste0(min(v), "-", max(v))
    } else {
      # If isolated (1) or just two (1,2), keep as comma separated
      paste(v, collapse = ",")
    }
  })
  
  # Combine groups with commas (e.g. "1-5, 8, 10-12")
  return(paste(ranges, collapse = ", "))
}

# 4. Generate the concise table
missing_months_concise <- CARIACO %>%
  select(date, all_of(vars)) %>%
  pivot_longer(cols = -date, names_to = "var_code", values_to = "value") %>%
  
  # Force explicit NAs for missing rows
  complete(var_code, date = full_dates) %>%
  
  mutate(
    Year = format(date, "%Y"),
    Month = month(date)
  ) %>%
  
  filter(Year >= 1995 & Year <= 2017) %>%
  
  group_by(var_code, Year) %>%
  
  # --- Apply custom range function ---
  summarize(
    missing_months = format_ranges(Month[is.na(value)]), 
    .groups = "drop"
  ) %>%
  
  # Join names and sort
  left_join(var_lookup, by = "var_code") %>%
  mutate(var_name = factor(var_name, levels = var_reorder)) %>%
  arrange(var_name) %>%
  
  # Final table shape
  select(var_name, Year, missing_months) %>%
  pivot_wider(names_from = Year, values_from = missing_months)

# 5. View and Export
print(missing_months_concise,n=21)


