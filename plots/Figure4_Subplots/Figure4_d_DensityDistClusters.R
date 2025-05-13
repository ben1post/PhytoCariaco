library(tidyverse)
require(cowplot)

library(cowplot)
library(ggpubr)

# get env data
CARIACO <- readRDS("data/processed/CARIACO_EnvData_combined.rds")


### ADD CLUSTER COMP ###

CARIACO$date <- as.Date(paste(CARIACO$time_month, "-15", sep=""), format="%m-%Y-%d")

# take negative of u10 (vector in u direction, mostly negative) as positive "wind speed":
CARIACO$u10_negative = -CARIACO$u10

extractGroups <- function(year1, year2, year3, year4){
  Group1 <- CARIACO %>%  
    filter(date >= as.Date(year1, format="%Y-%m-%d") & date <= as.Date(year2, format="%Y-%m-%d"))
  Group2 <- CARIACO %>%  
    filter(date >= as.Date(year3, format="%Y-%m-%d") & date <= as.Date(year4, format="%Y-%m-%d"))
  
  Group = rbind(Group1, Group2)
  
  return(Group)
}

# REMOVE DATA NOT IN CLUSTERING HERE

GroupRed = extractGroups(year1="1996-01-01", year2="2003-12-31", year3="2014-01-01", year4="2016-12-31") #2016-12-31
GroupBlue = extractGroups(year1="2004-01-01", year2="2013-12-31", year3="2017-06-01", year4="2017-12-31")


GroupRed$group = "Cluster 1 (1996-2003,2014-2016)"
GroupBlue$group = "Cluster 2 (2004-2013)"

ENV_DATA_groups = rbind(GroupRed, GroupBlue)

envdatmelt <- ENV_DATA_groups %>% 
  select(-date, -time_month) %>% 
  gather(variable, value, -group)

head(envdatmelt)



####### DENSITY DISTRIBUTION PLOTS ########

# create common legend:

colpallete = c("Cluster 1 (1996-2003,2014-2016)"="blue","Cluster 2 (2004-2013)"="red")

xxxx <- envdatmelt %>% filter(variable=="u10") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("10m U wind component") +theme_cowplot(font_size=20)+
  
  labs(color="Periods",fill="Periods") + theme(legend.position = "top") #+
  #plot_layout(tag_level = 'new')

legenddd <- get_legend(xxxx)
#plot_grid(legenddd)

A <- envdatmelt %>% filter(variable=="AMO") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("AMO") + theme(legend.position = "none") + xlab("[AMO]")

B <- envdatmelt %>% filter(variable=="MEIv2") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("MEI v2") + theme(legend.position = "none") + xlab("[MEI v2]")


C <- envdatmelt %>% filter(variable=="u10_negative") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete, legend="right") + ggtitle("Wind speed")+ theme(legend.position = "none") + xlab("[m/s]")

D <- envdatmelt %>% filter(variable=="tp") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("Precipitation") + theme(legend.position = "none") + xlab("[m Water/day]")

E <- envdatmelt %>% filter(variable=="e") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("Evaporation") + theme(legend.position = "none") + xlab("[m Water/day]")

prow1 = plot_grid(A, B, C,D,E, ncol=5,labels = c('d', 'e', 'f', 'g', 'h'))


F <- envdatmelt %>% filter(variable=="sst_10m") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("SST") + theme(legend.position = "none") + xlab("[°C]")

G <- envdatmelt %>% filter(variable=="Isotherm_21") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("21°C Isotherm Depth") + theme(legend.position = "none") + xlab("[m]")

H <- envdatmelt %>% filter(variable=="NO3_merged") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle(expression(NO[3])) + theme(legend.position = "none") + xlab("[µM]")

I <- envdatmelt %>% filter(variable=="PO4_merged") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete, legend="right") + ggtitle(expression(PO[4]))+ theme(legend.position = "none") + xlab("[µM]")

J <- envdatmelt %>% filter(variable=="SiO4_merged") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle(expression(SiO[4])) + theme(legend.position = "none") + xlab("[µM]")

prow2 = plot_grid(F, G, H, I, J, ncol=5, labels = c('i', 'j', 'k', 'l', 'm'))



K <- envdatmelt %>% filter(variable=="Salinity_bottles") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("Salinity") + theme(legend.position = "none") + xlab("[PSU]")

L <- envdatmelt %>% filter(variable=="Chlorophyll") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("Chlorophyll a") + theme(legend.position = "none") + xlab("[µM]") + scale_x_log10()


M <- envdatmelt %>% filter(variable=="GenusRichness") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete) + ggtitle("Genus Richness") + theme(legend.position = "none") + xlab("[Richness]")

N <- envdatmelt %>% filter(variable=="Shannon_gen") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete, legend="right") + ggtitle("Shannon Index")+ theme(legend.position = "none") + xlab("[Shannon Index]")

O <- envdatmelt %>% filter(variable=="Pielou_gen") %>%
  ggdensity(., x = "value",
            add = "median", rug = TRUE,
            color = "group", fill = "group",
            palette = colpallete, legend="right") + ggtitle("Pielou Index")+ theme(legend.position = "none") + xlab("[Pielou Index]")


prow3 = plot_grid(K, L, M, N, O, ncol=5,labels = c('n', 'o', 'p', 'q', 'r'))


options(repr.plot.width=12, repr.plot.height=8)
library(cowplot)

plottt <- plot_grid(legenddd, prow1,prow2,prow3, rel_heights = c(0.3,1,1,1), ncol=1, align='v', axis='r')

#ggsave("ComparativePlot2.pdf",plottt, width=12, height=11)





# CREATE STATISTICS TABLE

dataaa <- envdatmelt %>% filter(variable=="Isotherm_21")

wilcox.test(value ~ group, data = envdatmelt %>% filter(variable=="Isotherm_21"))

testout = wilcox.test(value ~ group, data = envdatmelt %>% filter(variable=="AMO"))

testout$p.value
testout$statistic

test_vars = c("AMO", "MEIv2", "u10", "tp", "e", "Isotherm_21",
              "sst_10m", "NO3_merged", "PO4_merged", "SiO4_merged",
              "Salinity_bottles", "Chlorophyll", "GenusRichness", "Shannon_gen", "Pielou_gen")

test_vars_finalnames = c("AMO", "MEI v.2", "Precipitation", "Evaporation", "u-component 10 m wind speed", "Isotherm 21 °C",
                         "SST", "NO3", "PO4", "SiO4",
                         "Salinity","Chlorophyll a",  "Genus Richness", "Shannon Index", "Pielou Index")

test_out_w = list()
test_out_p = list()
for (test_var in test_vars){
  test_out_w[test_var] = wilcox.test(value ~ group, data = envdatmelt %>% filter(variable==test_var))$statistic
  test_out_p[test_var] = wilcox.test(value ~ group, data = envdatmelt %>% filter(variable==test_var))$p.value
}

test_out_w$AMO
test_out_p$AMO

# statistics table!!

stars.pval <- function(x){
  stars <- c("***", "**", "*", "n.s.")
  var <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

library("kableExtra")
options(scipen=3)
statistics_table = map2_dfr(test_out_w, test_out_p, ~ tibble(W = .x, p.value = .y))
statistics_table$Vars <- test_vars
statistics_table$Variable <- test_vars_finalnames
statistics_table$Significance <- stars.pval(statistics_table$p.value)
statistics_table_formatted <- statistics_table %>% 
  mutate(p.value = signif(p.value, digits=3)) %>% 
  select(Variable, W, p.value, Significance)
kbl(statistics_table_formatted, "latex", booktabs = T, digits=3,
    caption="Wilcoxon sum rank test with continuity correction for variables between the two clusters. Values for in-situ data are interpolated across the top 100 meters for each individual monthly cruise sampling.")








