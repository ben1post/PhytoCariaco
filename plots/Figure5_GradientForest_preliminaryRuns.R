library(tidyverse)
library(vegan)
library(viridis)
library(cowplot)

## gradientForest package requires functional gfortran compiler to compile source (install "gcc" via homebrew on Mac)
# for author BPost it only worked on R 4.3.3 (Mac OS Sequoia, M1) with modified Makeconf file (Following comment by User "Shadow" here: https://stackoverflow.com/questions/76096681/macos-brew-system-r-packages-fail-to-install-with-emutls-w)
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)

# Read Phytoplankton Data:
phyto_counts <- readRDS("data/processed/PhytoplanktonInterpolatedCounts.RDS")
phyto_counts = phyto_counts[complete.cases(phyto_counts$counts),]

phyto_counts[phyto_counts$Genus=="Pseudo-nitzschia",]$Genus <- "Pseudo.nitzschia"

# get environmental data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


extractMatrixFix <- function(env_factors){
  ### Function to extract and match phytoplankton counts to environmental data for GF analysis
  ds_genus <- phyto_counts %>% 
    filter(TaxonRank == "Genus" | TaxonRank == "Species")%>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts))  %>%
    arrange(desc(date))
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  Mesh_genus$time_month =format(Mesh_genus$date, format="%m-%Y")
  new_Mesh_genus = Mesh_genus %>% select(-date)
  
  Env_dat <- CARIACO %>% 
    select("time_month", all_of(env_factors))
  
  
  CARIACO_dat_joined <- list(Env_dat, 
                             new_Mesh_genus) %>% 
    reduce(full_join, by = c("time_month")) %>% select(-time_month)
  
  
  Comp_Matrix <- CARIACO_dat_joined[complete.cases(CARIACO_dat_joined),]
  
  Matrix_Env <- Comp_Matrix %>% select(all_of(env_factors))
  Matrix_Genus <- Comp_Matrix %>% select(-all_of(env_factors))
  
  return(list(Matrix_Env, Matrix_Genus))
}


# Collect list of environmental factors
GF_env_factors_FULL = c("u10_lag1", "NO3_merged",
                        "PO4_merged", "SiO4_merged", "tp_lag1", "e_lag3",
                        "Salinity_bottles_lag3", "sst_10m", "Isotherm_21",
                        "AMO_lag2", 
                        "MEIv2_lag4")

GF_inputs <- extractMatrixFix(GF_env_factors_FULL)

envGF <- GF_inputs[[1]]
specGF <- GF_inputs[[2]]


mo5_specGF<- specGF %>% select_if(colSums(. > 0) > 10)

colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))

colnames(mo5_specGF)
colnames(envGF)

nSites <- dim(mo5_specGF)[1]
nSpecs <- dim(mo5_specGF)[2]
lev <- floor(log2(nSites * 0.368/2))
lev

print(nSites)
print(nSpecs)

dim(mo5_specGF)
dim(envGF)

# log transform species counts:
log_specGF <- log1p(mo5_specGF*100)

gf <- gradientForest(cbind(envGF, log_specGF),
                     predictor.vars=colnames(envGF),
                     response.vars=colnames(mo5_specGF),
                     ntree=500,
                     transform=NULL,
                     compact=T,
                     trace=T,
                     nbin=201,
                     maxLevel=lev,
                     corr.threshold=0.5)
gf


most_important <- names(importance(gf))[1:4]
par(mgp = c(2, 0.75,0))

# PLOT 1
plot(gf, plot.type="O")



#### START RUNS FOR TIME LAGS ####

GF_env_factors_NO3 = c("NO3_merged", "NO3_merged_lag1", "NO3_merged_lag2", "NO3_merged_lag3")

GF_env_factors_PO4 = c("PO4_merged", "PO4_merged_lag1", "PO4_merged_lag2", "PO4_merged_lag3")

GF_env_factors_SiO4 = c("SiO4_merged", "SiO4_merged_lag1", "SiO4_merged_lag2", "SiO4_merged_lag3")

GF_env_factors_SST = c("sst_10m", "sst_10m_lag1", "sst_10m_lag2", "sst_10m_lag3")

GF_env_factors_Isotherm_21 = c("Isotherm_21", "Isotherm_21_lag1", "Isotherm_21_lag2", "Isotherm_21_lag3")

GF_env_factors_Salinity_bottles = c("Salinity_bottles", "Salinity_bottles_lag1", "Salinity_bottles_lag2", "Salinity_bottles_lag3")

GF_env_factors_u10 = c("u10", "u10_lag1", "u10_lag2", "u10_lag3")

GF_env_factors_AMO = c("AMO", "AMO_lag1", "AMO_lag2", "AMO_lag3", "AMO_lag4", "AMO_lag4", "AMO_lag5", "AMO_lag6", "AMO_lag12")

GF_env_factors_MEIv2 = c("MEIv2", "MEIv2_lag1", "MEIv2_lag2", "MEIv2_lag3", "MEIv2_lag4", "MEIv2_lag4", "MEIv2_lag5", "MEIv2_lag6", "MEIv2_lag12")

GF_env_factors_tp = c("tp", "tp_lag1", "tp_lag2", "tp_lag3")

GF_env_factors_e = c("e", "e_lag1", "e_lag2", "e_lag3")


runGFmodel <- function(env_factors){
  # function to extract data and run GF model
  data <- extractMatrixFix(env_factors)
  
  envGF <- data[[1]]
  specGF <- data[[2]]
  
  mo5_specGF<- specGF %>% select_if(colSums(. > 0) > 10)
  colnames(mo5_specGF) <- make.names(colnames(mo5_specGF))
  
  nSites <- dim(mo5_specGF)[1]
  nSpecs <- dim(mo5_specGF)[2]
  lev <- floor(log2(nSites * 0.368/2))
  
  log_specGF <- log1p(mo5_specGF*100)
  
  gf <- gradientForest(cbind(envGF, log_specGF),
                       predictor.vars=colnames(envGF),
                       response.vars=colnames(mo5_specGF),
                       ntree=1500,
                       transform=NULL,
                       compact=T,
                       trace=T,
                       nbin=201,
                       maxLevel=lev,
                       corr.threshold=0.5)
  return(gf)
}



GF_output_NO3 <- runGFmodel(GF_env_factors_NO3)

GF_output_PO4 <- runGFmodel(GF_env_factors_PO4)

GF_output_SiO4 <- runGFmodel(GF_env_factors_SiO4)

GF_output_SST <- runGFmodel(GF_env_factors_SST)

GF_output_Isotherm_21 <- runGFmodel(GF_env_factors_Isotherm_21)

GF_output_Salinity_bottles <- runGFmodel(GF_env_factors_Salinity_bottles)

GF_output_u10 <- runGFmodel(GF_env_factors_u10)

GF_output_AMO <- runGFmodel(GF_env_factors_AMO)

GF_output_MEIv2 <- runGFmodel(GF_env_factors_MEIv2)

GF_output_tp <- runGFmodel(GF_env_factors_tp)

GF_output_e <- runGFmodel(GF_env_factors_e)


NO3_lags = as.numeric(importance(GF_output_NO3, sort=FALSE))
PO4_lags = as.numeric(importance(GF_output_PO4, sort=FALSE))
SiO4_lags = as.numeric(importance(GF_output_SiO4, sort=FALSE))
Temperature_lags = as.numeric(importance(GF_output_SST, sort=FALSE))
Isotherm_21_lags = as.numeric(importance(GF_output_Isotherm_21, sort=FALSE))
Salinity_bottles_lags = as.numeric(importance(GF_output_Salinity_bottles, sort=FALSE))
u10_lags = as.numeric(importance(GF_output_u10, sort=FALSE))
AMO_lags = as.numeric(importance(GF_output_AMO, sort=FALSE))
MEIv2_lags = as.numeric(importance(GF_output_MEIv2, sort=FALSE))
tp_lags = as.numeric(importance(GF_output_tp, sort=FALSE))
e_lags = as.numeric(importance(GF_output_e, sort=FALSE))

sq <- seq(8)
variable_selection = data.frame("lag" = c(0:6,12),"NO3"=NO3_lags[sq], 
                                "PO4"=PO4_lags[sq], 
                                "SiO4"=SiO4_lags[sq], 
                                "Temperature"=Temperature_lags[sq],               
                                "Isotherm_21"=Isotherm_21_lags[sq], 
                                "Salinity_bottles"=Salinity_bottles_lags[sq],
                                "Wind speed"=u10_lags[sq],
                                "Precipitation"=tp_lags[sq],
                                "Evaporation"=e_lags[sq], 
                                "AMO"=AMO_lags[sq],
                                "MEIv2"=MEIv2_lags[sq])


library("kableExtra")
options(scipen=3)
options(knitr.kable.NA = '-')

variable_selection %>%
  kbl("latex", booktabs = T, digits=3,
      caption="Individual Gradient Forest model runs for each variable and the corresponding time lags. 
        For in-situ variables, for which no coverage extends over the time series, 
        we tested all measurements up to a lag of 3 months. 
        For climate variables we tested up to a lag of 6 months. 
        Time lag with maximum importance per variable are highlighted in bold and chosen for the final model run.") #%>%

#cell_spec(2, bold = ifelse(. > 0.3, TRUE, FALSE)) 











### TEST: run EnvFactors individually to test significance ###


runEachEnvFactor <- function(listoffactors){
  output_imp = list()
  for (envfac in listoffactors){
    tryCatch({
      out_gf <- runGFmodel(envfac)
      output_imp[[envfac]] <- as.numeric(importance(out_gf, sort=FALSE))[[1]]},
      error = function(e) {
        # Handle the error
        cat("ERROR:", e$message, "\n")
        output_imp[[envfac]] <- NA
      }
    )
  }
  
  return(output_imp)
}

GF_output_NO3 <- runEachEnvFactor(GF_env_factors_NO3)

GF_output_PO4 <- runEachEnvFactor(GF_env_factors_PO4)

GF_output_SiO4 <- runEachEnvFactor(GF_env_factors_SiO4)

GF_output_SST <- runEachEnvFactor(GF_env_factors_SST)

GF_output_Isotherm_21 <- runEachEnvFactor(GF_env_factors_Isotherm_21)

GF_output_Salinity_bottles <- runEachEnvFactor(GF_env_factors_Salinity_bottles)

GF_output_u10 <- runEachEnvFactor(GF_env_factors_u10)

GF_output_AMO <- runEachEnvFactor(GF_env_factors_AMO)

GF_output_MEIv2 <- runEachEnvFactor(GF_env_factors_MEIv2)

GF_output_tp <- runEachEnvFactor(GF_env_factors_tp)

GF_output_e <- runEachEnvFactor(GF_env_factors_e)


NO3_lags = as.numeric(GF_output_NO3)
PO4_lags = as.numeric(GF_output_PO4)
SiO4_lags = as.numeric(GF_output_SiO4)
Temperature_lags = as.numeric(GF_output_SST)
Isotherm_21_lags = as.numeric(GF_output_Isotherm_21)
Salinity_bottles_lags = as.numeric(GF_output_Salinity_bottles)
u10_lags = as.numeric(GF_output_u10)
AMO_lags = as.numeric(GF_output_AMO)
MEIv2_lags = as.numeric(GF_output_MEIv2)
tp_lags = as.numeric(GF_output_tp)
e_lags = as.numeric(GF_output_e)

sq <- seq(8)
variable_selection = data.frame("lag" = c(0:6,12),"NO3"=NO3_lags[sq], 
                                "PO4"=PO4_lags[sq], 
                                "SiO4"=SiO4_lags[sq], 
                                "Temperature"=Temperature_lags[sq],               
                                "Isotherm_21"=Isotherm_21_lags[sq], 
                                "Salinity_bottles"=Salinity_bottles_lags[sq], 
                                "AMO"=AMO_lags[sq],
                                "MEIv2"=MEIv2_lags[sq],
                                "Wind speed"=u10_lags[sq],
                                "Precipitation"=tp_lags[sq],
                                "Evaporation"=e_lags[sq])