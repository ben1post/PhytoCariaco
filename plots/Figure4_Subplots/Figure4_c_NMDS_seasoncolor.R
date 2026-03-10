library(tidyverse)
require(cowplot)
library(vegan)
library(cowplot)
library(ggpubr)



# read interpolated phytoplankton cell count data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")

unique(phyto_counts$FuncGroup)


## prepare data for NMDS:
# get Genus Matrix:
ds_genus <- phyto_counts %>% 
  filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% #| TaxonRank == "Group") %>%
  group_by(Genus, date) %>%
  summarise(Total = sum(counts))  %>%
  arrange(desc(date))
tail(ds_genus, n=10)

ds_genus$year <- format(ds_genus$date, format="%Y")

ds_genus = ds_genus %>% select(-year)#%>% filter(year>=1996 & year<=2016) %>% select(-year)

Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)

Mesh_genus_noemtpyrows = Mesh_genus[rowSums(Mesh_genus[, -1])>0, ]


# get env data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


CARIACO$u10_negative = -CARIACO$u10

new_Mesh_genus = Mesh_genus_noemtpyrows
new_Mesh_genus$time_month =format(new_Mesh_genus$date, format="%m-%Y")


CARIACO_dat_joined <- list(CARIACO, 
                           new_Mesh_genus
) %>% 
  reduce(full_join, by = c("time_month"))

sel_env_factors = c("u10_negative","Salinity_bottles", "sst_10m", "Isotherm_21", "NO3_merged","PO4_merged","SiO4_merged", "AMO", "MEIv2", "tp")

firstSpec = as.numeric(which(names(CARIACO_dat_joined)=="Acanthoica"))
lastSpec = as.numeric(which(names(CARIACO_dat_joined)=="Zygosphaera"))
print(c(firstSpec,lastSpec))

Full_Matrix <- CARIACO_dat_joined %>%
  select("date", all_of(sel_env_factors), firstSpec:lastSpec)

Comp_Full_Matrix <- Full_Matrix[complete.cases(Full_Matrix),]


nfac = length(sel_env_factors)+2
nfac

Env_Matrix <- Comp_Full_Matrix[,sel_env_factors]

names(Env_Matrix)

names(Env_Matrix)[1] <- "Wind speed"
names(Env_Matrix)[2] <- "Salinity"

names(Env_Matrix)[3] <- "SST"

names(Env_Matrix)[4] <- "Isotherm Depth"


names(Env_Matrix)[5] <- "NO3"
names(Env_Matrix)[6] <- "PO4"
names(Env_Matrix)[7] <- "SiO4"


names(Env_Matrix)[8] <- "AMO"

names(Env_Matrix)[9] <- "MEIv2"


names(Comp_Full_Matrix)[nfac+222]
names(Comp_Full_Matrix)[nfac]

Genus_Matrix <- Comp_Full_Matrix[,nfac:(nfac+222)]
Genus_Matrix_2 <- Genus_Matrix %>% select_if(colSums(. > 0) > 5)
dim(Genus_Matrix)
dim(Genus_Matrix_2)

#convert com to a matrix
# add cube root transformation
m_com = as.matrix(Genus_Matrix_2^(1/3))

#nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray", k=2, trymax=200)
nmds

# ENVFIT
en = envfit(nmds, Env_Matrix, permutations = 999, na.rm = TRUE)

en

plot(nmds)
plot(en)


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores_2 = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores_2$Year = format(Comp_Full_Matrix$date, "%Y")

k_3 <- readRDS("plots/Figure4_Subplots/k_3.RDS")
data.scores_2$cluster = as.character(k_3[data.scores_2$Year])

head(data.scores_2)

options(repr.plot.width=10, repr.plot.height=10)

library(ggrepel)

#colsvec <- c("2" = "blue", "1" = "red")

# Add season column using your function
get_season_3cat <- function(dates) {
  months <- format(dates, "%m")
  season_lookup <- c(
    "12" = "Upwelling", "01" = "Upwelling", "02" = "Upwelling", 
    "03" = "Upwelling", "04" = "Upwelling", "05" = "Upwelling",
    "06" = "Secondary Upwelling", "07" = "Secondary Upwelling", "08" = "Secondary Upwelling",
    "09" = "Rainy", "10" = "Rainy", "11" = "Rainy"
  )
  factor(season_lookup[months], levels = c("Upwelling", "Secondary Upwelling", "Rainy"))
}

# Add season to data.scores_2
data.scores_2$season <- get_season_3cat(Comp_Full_Matrix$date)

# Define season colors
season_cols <- c("Upwelling" = "#7570b3", "Secondary Upwelling" = "#1b9e77", "Rainy" = "#d95f02")

# Plot with season coloring
NMDSplot_season <- ggplot(data.scores_2, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes(colour = season)) + 
  scale_colour_manual(values = season_cols) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(en_coord_cont), size = 4.5) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Season") +
  theme_cowplot(font_size = 20) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

NMDSplot_season


ggsave("plots/exports/FigureSX_seasonalNMDS_v1.pdf", NMDSplot_season, width=9, height=6)



###### PERMANOVA Analysis ########

# Get years for each row
years <- format(Comp_Full_Matrix$date, "%Y")

# Identify rows with valid cluster assignments (years present in k_3)
valid_rows <- years %in% names(k_3)

# Filter all objects to only valid rows
m_com_filtered <- m_com[valid_rows, ]
cluster_factor <- as.factor(k_3[years[valid_rows]])
season_factor <- get_season_3cat(Comp_Full_Matrix$date[valid_rows])

# Check how many rows removed
cat("Rows removed:", sum(!valid_rows), "\n")
cat("Rows remaining:", sum(valid_rows), "\n")
# PERMANOVA - testing both term orders to check for collinearity
# Order 1: Cluster first
set.seed(123)
perm_clust_first <- adonis2(m_com_filtered ~ cluster_factor + season_factor, 
                            method = "bray", 
                            permutations = 999,
                            by = "terms")

# Order 2: Season first  
set.seed(123)
perm_season_first <- adonis2(m_com_filtered ~ season_factor + cluster_factor, 
                             method = "bray", 
                             permutations = 999,
                             by = "terms")

# View results
print("=== Cluster entered first ===")
print(perm_clust_first)

print("=== Season entered first ===")
print(perm_season_first)



# =============================================================================
# LATEX TABLE OUTPUT FOR PERMANOVA RESULTS
# =============================================================================

format_p_latex <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("$<$0.001")
  return(sprintf("%.3f", p))
}

format_r2 <- function(r2) {
  sprintf("%.1f", r2 * 100)
}

format_f <- function(f) {
  if (is.na(f)) return("—")
  sprintf("%.1f", f)
}

# Extract values from the PERMANOVA output (cluster entered first)
perm_df <- as.data.frame(perm_clust_first)

# Build table as a single string
latex_table <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Permutational multivariate analysis of variance (PERMANOVA) testing the effects of ",
  "interannual cluster and season on phytoplankton community composition. ",
  "Community data (genus-level cell counts) were cube-root transformed prior to analysis. ",
  "Bray--Curtis dissimilarity was used as the distance metric, with 999 permutations. ",
  "Observations from 1995 and 2017 were excluded due to missing cluster assignments ($n$ = 193). ",
  "Results shown are for cluster entered first; reversing term order yielded nearly identical ",
  "R\\textsuperscript{2} values (cluster: 17.9\\%, season: 3.5\\%), indicating minimal collinearity between factors.}\n",
  "\\label{sup:tab:permanova_season_cluster}\n",
  "\\small\n",
  "\\begin{tabular}{l rrrrl}\n",
  "\\hline\n",
  "Factor & df & Sum of Squares & R\\textsuperscript{2} (\\%) & \\textit{F} & \\textit{p}-value \\\\\n",
  "\\hline\n"
)

# Add data rows
row_names <- c("Cluster", "Season", "Residual", "Total")
for (i in 1:4) {
  if (i <= 2) {
    # Cluster and Season rows
    latex_table <- paste0(latex_table,
                          sprintf("%s & %d & %.3f & %s & %s & %s \\\\\n",
                                  row_names[i],
                                  perm_df$Df[i],
                                  perm_df$SumOfSqs[i],
                                  format_r2(perm_df$R2[i]),
                                  format_f(perm_df$F[i]),
                                  format_p_latex(perm_df$`Pr(>F)`[i])))
  } else if (i == 3) {
    # Residual row (no F or p-value)
    latex_table <- paste0(latex_table,
                          sprintf("%s & %d & %.3f & %s & — & — \\\\\n",
                                  row_names[i],
                                  perm_df$Df[i],
                                  perm_df$SumOfSqs[i],
                                  format_r2(perm_df$R2[i])))
  } else {
    # Total row (no F or p-value)
    latex_table <- paste0(latex_table,
                          sprintf("%s & %d & %.3f & %s & — & — \\\\\n",
                                  row_names[i],
                                  perm_df$Df[i],
                                  perm_df$SumOfSqs[i],
                                  format_r2(perm_df$R2[i])))
  }
}

# Close table
latex_table <- paste0(latex_table,
                      "\\hline\n",
                      "\\end{tabular}\n",
                      "\\end{table}\n"
)

# Output complete table
cat("\n=== LATEX TABLE ===\n\n")
cat(latex_table)

