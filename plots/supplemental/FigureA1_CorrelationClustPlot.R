source("plots/Figure2_ZScoreTimeSeries.R")

#### Create large corrplot of all variables for visualisation and supplemental! ####

library(corrplot)
envvars <- CARIACO %>% select(MEIv2, AMO, u10_negative, tp, e, 
                              sst_10m, Isotherm_21, Salinity_bottles, NO3_merged, PO4_merged, SiO4_merged,
                              Chlorophyll, Abundance_Diatom, Abundance_Hapto, Abundance_Dino, Abundance_Cyano, Abundance_Nanoflagellate,
                              GenusRichness, Shannon_gen, Pielou_gen,
                              -time_month)

names(envvars) <- c("MEI v.2", "AMO", "Wind speed", "Precipitation", "Evaporation", 
                    "SST", "21°C Isotherm", "Salinity", "NO3", "PO4", "SiO4",
                    "Chlorophyll a", "Diatoms", "Haptophytes", 
                    "Dinoflagellates", "Cyanobacteria", "Nanoflagellates",
                    "Genus Richness", "Shannon Index", "Pielou's Index")

m = cor(envvars, use = "pairwise.complete.obs", method="spearman")
#options(repr.plot.width=15, repr.plot.height=15)
corrplot(m, diag = FALSE, order = 'hclust', addrect = 2)

testRes = cor.mtest(envvars, conf.level = 0.95, method="spearman")
corrplot(m, p.mat = testRes$p, sig.level = 0.10, order = 'hclust', addrect = 4, insig='blank')


pdf(file = "plots/exports/CorrClustSupplementalPlot_v1.pdf", width = 6, height=6)

corrplot(m,diag = FALSE, p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', tl.col = 'black')

corrplot(m,diag = FALSE,order = 'hclust', addrect = 5, p.mat = testRes$p, sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', tl.col = 'black')

dev.off()
