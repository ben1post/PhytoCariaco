library(tidyr)
library(dplyr)
library(lubridate)

# Read data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")
CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")

# Transform variables
CARIACO$u10_negative <- -CARIACO$u10
CARIACO$Chlorophyll_log <- log1p(CARIACO$Chlorophyll)
CARIACO$PrimaryProductivity_log <- log1p(CARIACO$PrimaryProductivity)
CARIACO$Abundance_Diatom_log <- log1p(CARIACO$Abundance_Diatom)
CARIACO$Abundance_Hapto_log <- log1p(CARIACO$Abundance_Hapto)
CARIACO$Abundance_Dino_log <- log1p(CARIACO$Abundance_Dino)
CARIACO$Abundance_Cyano_log <- log1p(CARIACO$Abundance_Cyano)
CARIACO$Abundance_Nanoflagellate_log <- log1p(CARIACO$Abundance_Nanoflagellate)

# Variables
vars <- c("MEIv2", "AMO", "u10_negative", "tp", "Isotherm_21", "Salinity_bottles", "sst_10m",
          "NO3_merged", "PO4_merged", "SiO4_merged", "Chlorophyll_log", "PrimaryProductivity_log",
          "Abundance_Diatom_log", "Abundance_Hapto_log", "Abundance_Dino_log", 
          "Abundance_Cyano_log", "Abundance_Nanoflagellate_log",
          "GenusRichness", "Shannon_gen", "Pielou_gen")

renamed_vars <- c("MEIv2", "AMO", "Wind speed", "Precipitation", "Isotherm Depth", "Salinity", "SST",
                  "NO3", "PO4", "SiO4", "Chlorophyll a", "Primary Productivity",
                  "Diatoms", "Haptophytes", "Dinoflagellates", "Cyanobacteria", "Nanoflagellates",
                  "Richness", "Shannon", "Pielou")

var_reorder <- c("MEIv2", "AMO", "Wind speed", "Precipitation", "SST", "Isotherm Depth", "Salinity",
                 "NO3", "PO4", "SiO4", "Chlorophyll a", "Primary Productivity", "Diatoms", 
                 "Haptophytes", "Dinoflagellates", "Cyanobacteria", "Nanoflagellates",
                 "Richness", "Shannon", "Pielou")

var_lookup <- data.frame(var_code = vars, var_name = renamed_vars, stringsAsFactors = FALSE)

full_dates <- seq(as.Date("1995-01-15"), as.Date("2017-12-15"), by = "month")

# Format ranges function
format_ranges <- function(x) {
  if (length(x) == 0) return("-")
  x <- sort(unique(x))
  groups <- cumsum(c(1, diff(x) != 1))
  split_x <- split(x, groups)
  ranges <- sapply(split_x, function(v) {
    if (length(v) > 2) paste0(min(v), "-", max(v))
    else paste(v, collapse = ",")
  })
  paste(ranges, collapse = ",")
}

# Generate table
missing_table <- CARIACO %>%
  select(date, all_of(vars)) %>%
  pivot_longer(cols = -date, names_to = "var_code", values_to = "value") %>%
  complete(var_code, date = full_dates) %>%
  mutate(Year = format(date, "%Y"), Month = month(date)) %>%
  filter(Year >= 1995 & Year <= 2017) %>%
  group_by(var_code, Year) %>%
  summarize(missing = format_ranges(Month[is.na(value)]), .groups = "drop") %>%
  left_join(var_lookup, by = "var_code") %>%
  mutate(var_name = factor(var_name, levels = var_reorder)) %>%
  arrange(var_name) %>%
  select(var_name, Year, missing) %>%
  pivot_wider(names_from = Year, values_from = missing)

print(missing_table, n = 21, width = Inf)

# LaTeX output
library(kableExtra)

latex_output <- missing_table %>%
  rename(Variable = var_name) %>%
  kbl(format = "latex", booktabs = TRUE, linesep = "",
      caption = "Missing months per variable and year (1 = January, 12 = December; - = complete).",
      label = "tab:missing-data") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"), font_size = 8) %>%
  row_spec(0, angle = 90)


# LaTeX output - landscape table, with collapsed microscopy rows

# Collapse microscopy variables into one row
microscopy_vars <- c("Diatoms", "Haptophytes", "Dinoflagellates", "Cyanobacteria", 
                     "Nanoflagellates", "Richness", "Shannon", "Pielou")

missing_table_collapsed <- missing_table %>%
  mutate(var_name = as.character(var_name)) %>%
  mutate(var_name = ifelse(var_name %in% microscopy_vars, "Microscopy cell counts", var_name)) %>%
  distinct()

# Generate LaTeX
years <- names(missing_table_collapsed)[-1]
n_years <- length(years)

header_row <- paste0(paste(years, collapse = " & "), " \\\\")

data_rows <- apply(missing_table_collapsed, 1, function(row) {
  paste0(row[1], " & ", paste(row[-1], collapse = " & "), " \\\\")
})

latex_table <- c(
  "\\begin{sidewaystable}",
  "\\centering",
  "\\caption{Missing months per variable and year in the CARIACO time series (1995-2017). Numbers indicate missing months (1 = January, 12 = December); ranges indicate consecutive missing months (e.g., 1-10 = January through October); `-' indicates complete data for that year.}",
  "\\label{tab:missing-data}",
  "\\resizebox{\\textheight}{!}{%",
  paste0("\\begin{tabular}{l", paste(rep("c", n_years), collapse = ""), "}"),
  "\\hline",
  paste0("Variable & ", header_row),
  "\\hline",
  data_rows,
  "\\hline",
  "\\end{tabular}}",
  "\\end{sidewaystable}"
)

cat(paste(latex_table, collapse = "\n"))

