# Script to merge all processed datasets into cohesive data frame for analysis
# All lag calculations are performed here after creating a complete monthly sequence

library(tidyverse)

options(stringsAsFactors = FALSE)

# ============================================
# READ ALL DATASETS (without lag calculations)
# ============================================

# ERA5 Data
wind_ds <- readRDS("data/processed/ERA5_data.rds") %>%
  mutate(
    date = as.Date(time),
    time_month = format(date, format = "%m-%Y")
  ) %>%
  select(time_month, u10, e, tp)

# Niskin Data (NOTE: Chl, PP, nutrients are depth-integrated)
niskin_ds <- readRDS("data/processed/Niskin_qchecked_100m.RDS") %>%
  mutate(time_month = format(date, format = "%m-%Y")) %>%
  select(-date)

niskin_ds$PrimaryProductivity <- niskin_ds$PrimaryProductivity * 12

# Phytoplankton Data (NOTE: abundances are depth-integrated)
diversity_ds <- readRDS("data/processed/PhytoplanktonAbundanceDiversity.RDS")

# CTD Data
ctd_ds <- readRDS("data/processed/CTD_Isotherm21_MLD.rds")

# AMO Data
amo_table <- read.table("data/AMO/amo_monthly.txt", header = FALSE, skip = 1)
nm2 <- setNames(sprintf("%02d", 1:12), paste0("V", 2:13))

amo_ds <- amo_table %>%
  gather(Month, AMO, -V1) %>%
  mutate(
    Month = nm2[Month],
    time_month = paste(Month, V1, sep = "-")
  ) %>%
  select(time_month, AMO)

# MEI v.2 Data
meiv2_lines <- readLines("data/MEIv2/meiv2.data")
meiv2_lines <- meiv2_lines[c(-1, -48:-51)]

meiv2_ds <- read.table(textConnection(meiv2_lines), header = FALSE) %>%
  gather(Month, MEIv2, -V1) %>%
  mutate(
    Month = nm2[Month],
    time_month = paste(Month, V1, sep = "-")
  ) %>%
  select(time_month, MEIv2)

# ============================================
# CREATE COMPLETE MONTHLY BACKBONE
# ============================================

# Define CARIACO time series period
start_date <- as.Date("1995-11-01")
end_date <- as.Date("2017-01-01")

# Extended backbone (to accommodate lag variables of climate indices)
backbone_start <- as.Date("1985-11-01")  # 10 years before study start

extended_months <- data.frame(
  date = seq(from = backbone_start, to = end_date, by = "month")
) %>%
  mutate(time_month = format(date, format = "%m-%Y"))

cat("Extended backbone:", nrow(extended_months), "months\n")
cat("  From:", format(backbone_start, "%b %Y"), "to", format(end_date, "%b %Y"), "\n")


# Create complete monthly sequence
complete_months <- data.frame(
  date = seq(from = start_date, to = end_date, by = "month")
) %>%
  mutate(time_month = format(date, format = "%m-%Y"))

cat("Complete monthly backbone:", nrow(complete_months), "months\n")

# ============================================
# JOIN ALL DATASETS TO MONTHLY BACKBONE
# ============================================

CARIACO_extended <- extended_months %>%
  left_join(wind_ds, by = "time_month") %>%
  left_join(niskin_ds, by = "time_month") %>%
  left_join(diversity_ds, by = "time_month") %>%
  left_join(ctd_ds, by = "time_month") %>%
  left_join(amo_ds, by = "time_month") %>%
  left_join(meiv2_ds, by = "time_month") %>%
  arrange(date)

# ============================================
# CALCULATE ALL LAGS ON COMPLETE SEQUENCE
# ============================================

CARIACO_extended_lags <- CARIACO_extended %>%
  mutate(
    # ERA5 wind/climate lags (1-6 months)
    across(
      c(u10, e, tp),
      list(lag1 = ~lag(., 1), lag2 = ~lag(., 2), lag3 = ~lag(., 3),
           lag4 = ~lag(., 4), lag5 = ~lag(., 5), lag6 = ~lag(., 6))
    ),
    # Nutrient lags (1-3 months)
    across(
      c(NO3_merged, PO4_merged, SiO4_merged),
      list(lag1 = ~lag(., 1), lag2 = ~lag(., 2), lag3 = ~lag(., 3))
    ),
    # Physical oceanography lags (1-3 months)
    across(
      c(Temperature, Salinity_bottles),
      list(lag1 = ~lag(., 1), lag2 = ~lag(., 2), lag3 = ~lag(., 3))
    ),
    # CTD-derived lags (1-6 months)
    across(
      c(Isotherm_21, sst_10m),
      list(lag1 = ~lag(., 1), lag2 = ~lag(., 2), lag3 = ~lag(., 3),
           lag4 = ~lag(., 4), lag5 = ~lag(., 5), lag6 = ~lag(., 6))
    ),
    # Climate indices lags (1-6 months + annual)
    across(
      c(AMO, MEIv2),
      list(lag1 = ~lag(., 1), lag2 = ~lag(., 2), lag3 = ~lag(., 3),
           lag4 = ~lag(., 4), lag5 = ~lag(., 5), lag6 = ~lag(., 6),
           lag7 = ~lag(., 7), lag8 = ~lag(., 8), lag9 = ~lag(., 9),
           lag10 = ~lag(., 10), lag11 = ~lag(., 11), lag12 = ~lag(., 12),
           # --- YEAR LAGS ADDED HERE ---
           lag24 = ~lag(., 24),  # 2 Years
           lag36 = ~lag(., 36),  # 3 Years
           lag48 = ~lag(., 48),  # 4 Years
           lag60 = ~lag(., 60),   # 5 Years
           lag72 = ~lag(., 72),   # 6 Years
           lag84 = ~lag(., 84),   # 7 Years
           lag96 = ~lag(., 96),   # 8 Years
           lag108 = ~lag(., 108)   # 9 Years
          )
    ),
  )

# ============================================
# TRUNCATE TO STUDY PERIOD
# ============================================

CARIACO_dat_joined <- CARIACO_extended_lags %>%
  filter(date >= start_date & date <= end_date) %>%
  select(-date)

# ============================================
# DIAGNOSTICS
# ============================================

cat("\nFinal dataset:", nrow(CARIACO_dat_joined), "rows\n\n")

# Check NA patterns for key variables (excluding lag columns)
key_vars <- c("u10", "tp", "Chlorophyll", "PrimaryProductivity", 
              "NO3_merged", "PO4_merged", "SiO4_merged",
              "Isotherm_21", "sst_10m", "MLD",
              "Shannon_gen", "GenusRichness",
              "AMO", "MEIv2")

na_summary <- CARIACO_dat_joined %>%
  select(any_of(key_vars)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  mutate(pct_NA = round(100 * n_NA / nrow(CARIACO_dat_joined), 1)) %>%
  arrange(desc(n_NA))

cat("Missing data summary (key variables):\n")
print(na_summary)

# ============================================
# SAVE
# ============================================

saveRDS(CARIACO_dat_joined, "data/processed/CARIACO_EnvData_combined.rds")

cat("\nData saved to: data/processed/CARIACO_EnvData_combined.rds\n")
