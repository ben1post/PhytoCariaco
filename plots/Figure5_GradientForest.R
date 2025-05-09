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
    filter(TaxonRank == "Genus" | TaxonRank == "Species") %>% 
    group_by(Genus, date) %>%
    summarise(Total = sum(counts))  %>%
    arrange(desc(date))
  
  
  Mesh_genus <- pivot_wider(ds_genus, names_from = Genus, values_from = Total, values_fn = sum, values_fill = 0)
  
  new_Mesh_genus = Mesh_genus
  new_Mesh_genus$time_month =format(new_Mesh_genus$date, format="%m-%Y")
  
  CARIACO_dat_joined <- list(CARIACO, 
                             new_Mesh_genus
  ) %>% 
    reduce(full_join, by = c("time_month"))
  
  Matrix_02 <- CARIACO_dat_joined
  
  Matrix_03 <- Matrix_02 %>%
    select("date", all_of(env_factors), length(names(CARIACO)):length(Matrix_02))
  
  Comp_Matrix_04 <- Matrix_03[complete.cases(Matrix_03),]

  nfac = length(env_factors)+3
  
  Matrix_05_Env <- Comp_Matrix_04[,env_factors]
  Matrix_05_Genus <- Comp_Matrix_04[,nfac:(nfac+length(Mesh_genus[-1])-1)]
  
  return(list(Matrix_05_Env, Matrix_05_Genus))
}


# Collect list of environmental factors
GF_env_factors_FULL = c("u10", "Isotherm_21", "NO3_merged",
                        "PO4_merged", "SiO4_merged",
                        "Salinity_bottles", "Temperature",
                        "AMO_lag5", 
                        "MEIv2_lag4")

GF_inputs <- extractMatrixFix(GF_env_factors_FULL)

envGF <- GF_inputs[[1]]
specGF <- GF_inputs[[2]]


mo5_specGF<- specGF %>% select_if(colSums(. > 0) > 5)

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
                     ntree=1500,
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



# PLOT 2
options(repr.plot.width=7, repr.plot.height=5)

#pdf(file = "TestGFPlot_1.pdf",   # The directory you want to save the file in
#    width = 10, # The width of the plot in inches
#    height = 10) # The height of the plot in inches

grid.arrange(plot(gf, plot.type="O"),
             plot(gf, plot.type="C", imp.vars=most_important,
                  show.species = F, common.scale=T, cex.axis=0.6, cex.lab=0.7,
                  line.ylab=0.9, par.args = list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0))),
             ncol=2, nrow=1
)

#dev.off()


# PLOT 3
options(repr.plot.width=10, repr.plot.height=7)

#pdf(file = "TestGFPlot.pdf",   # The directory you want to save the file in
#    width = 10, # The width of the plot in inches
#    height = 7) # The height of the plot in inches


plot(gf, plot.type="S", imp.vars=most_important,
     leg.posn="topright", cex.legend=0.7, cex.axis=0.9,
     cex.lab=1.1, line.ylab=0.9, par.args=list(mgp=c(1.5,0.5,0), mar=c(3.1,1.5,0.1,1)))


#dev.off()


# PLOT 4
options(repr.plot.width=10, repr.plot.height=13)

#pdf(file = "TestGFPlot_3.pdf",   # The directory you want to save the file in
#    width = 10, # The width of the plot in inches
#    height = 7) # The height of the plot in inches

plot(gf, plot.type="C", imp.vars=most_important,
     show.overall = F, legend=T, leg.posn="topleft",
     leg.nspecies=5, cex.lab=0.7, cex.legend=0.7, 
     cex.axis=1,line.ylab=0.9,
     par.args = list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))

#dev.off()



# PLOT 5
#pdf(file = "TestGFPlot_2.pdf",   # The directory you want to save the file in
#    width = 5, # The width of the plot in inches
#    height = 10) # The height of the plot in inches

plot(gf, plot.type="C", imp.vars=most_important,
     show.species = F, common.scale=T, cex.axis=0.6, cex.lab=0.7,
     line.ylab=0.9, par.args = list(mgp=c(1.5,0.5,0), mar=c(2.5,1,0.1,0.5), omi=c(0,0.3,0,0)))

#dev.off()




###### FINAL PLOT OUTPUT ######

show.species=FALSE
show.overall=TRUE

imp.vars <- imp.var.names <- names(importance(gf))[1:9]
imp.vars

par(mfrow=rev(n2mfrow(length(imp.vars)*(show.species+show.overall))))
cols <- rainbow(length(names(gf$result)))
names(cols) <- names(gf$result)

xaxt <- if(show.overall) "n" else "s"



las = 1
cex.axis = 0.7
cex.names = cex.axis
horiz = TRUE

#imp.a <- importance(gf,"Accuracy")
imp.w <- importance(gf,"Weighted")
#o.a <- order(imp.a)
o.w <- order(imp.w)
varnames = names(imp.w[o.w])
data = imp.w[o.w]
data_df <- data.frame(var=factor(names(imp.w), levels=names(imp.w)), importance=imp.w, row.names=NULL)
title= expression(paste(R^2, " weighted importance"))

weightedImp_plot <- ggplot() + geom_bar(data=data_df[o.w,], aes(x=var, y=importance), stat="identity") +theme_cowplot()+
  coord_flip() + scale_x_discrete(limits = rev) + xlab("Predictor Variable") + ylab(title)

weightedImp_plot


imp.vars.names=imp.vars
common.scale=F
line.ylab=1.0
cex.legend=0.75

#species
leg.nspecies=10
leg.posn="topleft"
legend=TRUE


linesdat = list()

for (varX in imp.vars) {
  print(varX)
  CU <- cumimp(gf, varX, "Species")
  xlim <- range(sapply(CU, "[[", "x"))
  ylim <- range(sapply(CU, "[[", "y"))
  
  specdat = list()
  
  isub <- seq(1,length(CU[[species]]$x),len=pmin(500,length(CU[[species]]$x)))
  
  xvals = CU[[species]]$x[isub]
  yvals = CU[[species]]$y[isub]
  
  columns = c("Predictor","Species","x","y")
  testdf = data.frame(matrix(nrow=length(xvals), ncol=length(columns)))
  names(testdf) <- columns
  
  testdf$x <- xvals
  testdf$y <- yvals
  testdf$Species <- species
  testdf$Predictor <- varX
  
  specdat[[species]] = testdf
  
  linesdat[[varX]] = rbindlist(specdat)
  
  no.species<-length(names(cols))
  # only label most important species
  imp.sp <- sapply(CU, function(cu) max(cu$y))
  best <- order(-imp.sp)[1:min(leg.nspecies,length(imp.sp))]
}

merged_spec_cumimp_df <- rbindlist(linesdat)


length(unique(merged_spec_cumimp_df$Species))


specs = unique(merged_spec_cumimp_df$Species)
specs

# func group color scale:

GenusFuncGroup_match <- phyto_counts %>% group_by(Genus) %>% summarize(FuncGroup = first(FuncGroup))

GenusMatched <- GenusFuncGroup_match[GenusFuncGroup_match$Genus %in% specs,]

GenusMatched

# calculate the proportion of counts that are predicted by GF:
sum(specGF[GenusMatched$Genus])/ sum(specGF) # captures this proportion of counts!


library(Polychrome)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

genscolscale = gg_color_hue(length(unique(GenusMatched$FuncGroup)))
gens = unique(GenusMatched$FuncGroup)
names(genscolscale) <- gens
genscolscale



options(repr.plot.width=2, repr.plot.height=2)

TestGenScalePlot <- ggplot(data = data.frame(FuncGroup = names(genscolscale), Color=genscolscale, dat=1:5, row.names = NULL)) + 
  geom_point(aes(x=FuncGroup, y=dat, col=FuncGroup)) + theme_cowplot() + scale_colour_manual(values = genscolscale)

ColScaleLegend <- get_legend(TestGenScalePlot)
plot_grid(ColScaleLegend)

GenusMatched$color <-  genscolscale[GenusMatched$FuncGroup]


# COLOR SCALE:
options(repr.plot.width=6, repr.plot.height=5)
speccolscale = as.character(GenusMatched$color)
names(speccolscale) <- as.character(GenusMatched$Genus)

swatch(speccolscale)


# SPEC CUM IMP PLOT
str(merged_spec_cumimp_df)
merged_spec_cumimp_df$Predictor = factor(merged_spec_cumimp_df$Predictor, levels=names(imp.w))

specCumImp_Plot <- ggplot(data=merged_spec_cumimp_df)+geom_line(aes(x=x,y=y, col=Species), linetype = "dashed")+
  facet_wrap(~Predictor,scales = "free_x") + xlab("") + ylab("") + theme_cowplot(font_size=12) + scale_colour_manual(values = speccolscale) + guides(color="none")
specCumImp_Plot

# OVERALL CUM IMP PLOT
# Overall Cum Imp Plot
#
# Combine species
#
plotlist = list()
overallcumimp = list()

for (varX in imp.vars) {
  print(varX)
  CU <- cumimp(gf, varX)
  ymax <- max(CU$y)
  if (varX == imp.vars[1]){ymax1 <- ymax}
  isub <- seq(1,length(CU$x),len=pmin(500,length(CU$x)))

  xvals=CU$x[isub]
  yvals=CU$y[isub]
  
  columns = c("Predictor","Species","x","y")
  testdf = data.frame(matrix(nrow=length(xvals), ncol=length(columns)))
  names(testdf) <- columns
  
  testdf$x <- CU$x[isub]
  testdf$y <- CU$y[isub]
  testdf$Species <- "overall"
  testdf$Predictor <- varX
  
  overallcumimp[[varX]] = testdf
}

OVCumImpDat <- rbindlist(overallcumimp)


OVCumImpDat$Predictor = factor(OVCumImpDat$Predictor, levels=names(imp.w))


# GET GENUS LEVEL PREDICTORS
linesdat_genus = list()

for (varX in imp.vars) {
  print(varX)
  CU <- cumimp(gf, varX, "Species")
  xlim <- range(sapply(CU, "[[", "x"))
  ylim <- range(sapply(CU, "[[", "y"))
  
  importance.df <- gf$res[gf$res$var==varX,]
  
  FGdat = list()
  for(funcgroup in unique(GenusMatched$FuncGroup)) {
    
    Genus <- GenusMatched %>% filter(FuncGroup==funcgroup) %>% pull(Genus)
    
    
    Gen_CU <- getCU(subset(importance.df, spec %in% Genus | var == varX), mean(x$imp.rsq[varX, Genus]))
    
    isub <- seq(1,length(Gen_CU$x),len=pmin(500,length(Gen_CU$x)))
    
    xvals = Gen_CU$x[isub]
    yvals = Gen_CU$y[isub]
    
    columns = c("Predictor","FuncGroup","x","y")
    testdf = data.frame(matrix(nrow=length(xvals), ncol=length(columns)))
    names(testdf) <- columns
    
    testdf$x <- xvals
    testdf$y <- yvals
    testdf$FuncGroup <- funcgroup
    testdf$Predictor <- varX
    
    FGdat[[funcgroup]] = testdf
  }
  
  linesdat_genus[[varX]] = rbindlist(FGdat)
}

merged_genus_cumimp_df <- rbindlist(linesdat_genus)


str(merged_genus_cumimp_df$Predictor)
str(OVCumImpDat$Predictor)

merged_genus_cumimp_df$Predictor <- factor(merged_genus_cumimp_df$Predictor, levels=levels(OVCumImpDat$Predictor))

options(repr.plot.width=10, repr.plot.height=12)

overallCumImp_Plot <- ggplot(data=merged_genus_cumimp_df)+ geom_line(aes(x=x,y=y, col=FuncGroup), linetype = "dashed", alpha=1)+
  geom_line(data=OVCumImpDat, aes(x=x,y=y))+
  facet_wrap(~Predictor,scales = "free_x") + xlab("") + ylab("") + theme_cowplot(font_size=16) + scale_colour_manual(values = genscolscale) + guides(color="none")
overallCumImp_Plot


#### SPECIES IMP PLOT ####

options(repr.plot.width=7, repr.plot.height=7)

perf <- importance(gf, type="Species") 
#perf

o.s <- order(perf)
specnames = names(perf[o.s])
speccdata = perf[o.s]

speccdata_df <- data.frame(Genus=factor(specnames, levels=specnames), importance=speccdata, row.names=NULL)

n <- length(perf)
title= expression(paste(R^2, " weighted importance"))

specImp_plot <- ggplot(data=speccdata_df) + geom_point(aes(x=Genus,y=importance, col=Genus)) + theme_cowplot(font_size=10) + 
  coord_flip() + ylab(title) + scale_colour_manual(values = speccolscale) + guides(color="none")


specImp_plot 



##### COMBINED PLOT #####

options(repr.plot.width=12, repr.plot.height=16)

left_column = overallCumImp_Plot
right_column = plot_grid(weightedImp_plot, specImp_plot,ColScaleLegend, ncol=3, rel_widths=c(1,1,0.1))

GF_output_plot1 <- plot_grid(right_column, left_column, ncol=1, rel_heights = c(1.8,3))

GF_output_plot1
#ggsave("GF_output_plot3.pdf",GF_output_plot1, width=12, height=16)


