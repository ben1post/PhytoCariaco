library(tidyverse)
require(cowplot)

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

ds_genus = ds_genus %>% filter(year>=1996 & year<=2016) %>% select(-year)

Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)

Mesh_genus_noemtpyrows = Mesh_genus[rowSums(Mesh_genus[, -1])>0, ]


# get env data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


new_Mesh_genus = Mesh_genus_noemtpyrows
new_Mesh_genus$time_month =format(new_Mesh_genus$date, format="%m-%Y")


CARIACO_dat_joined <- list(CARIACO, 
                           new_Mesh_genus
) %>% 
  reduce(full_join, by = c("time_month"))

sel_env_factors = c("u10","Salinity_bottles", "Temperature", "Isotherm_21", "NO3_merged","PO4_merged","SiO4_merged", "AMO_lag2", "MEIv2_lag4")#, "Shannon_gen", "Pielou_gen", "GenusRichness")

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

names(Env_Matrix)[1] <- "u-comp wind 10m"
names(Env_Matrix)[2] <- "Salinity"

names(Env_Matrix)[3] <- "Temperature"

names(Env_Matrix)[4] <- "Isotherm_21"


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

colsvec <- c("2" = "blue", "1" = "red")

NMDSplot <- ggplot(data.scores_2, aes(x = NMDS1, y = NMDS2)) + 
  
  
  geom_point(size = 2, aes(colour = cluster))+ 
  
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  
  labs(x = "NMDS1", colour = "Cluster", y = "NMDS2", shape = "Group")  + 
  
  scale_colour_manual(values = colsvec) + 
  
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth =1, alpha = 0.5, colour = "grey30") +
  
  
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
                  fontface = "bold", label = row.names(en_coord_cont), size=4.5) + theme_cowplot((font_size=20)) + guides(colour="none")

NMDSplot

