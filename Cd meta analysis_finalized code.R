library(ggplot2)
library(dplyr)
library(maps)
library(ggpubr)
library(raster)
library(sp)
library(readr)
library(gridExtra)
library(ggrepel)
library(eulerr)
library(gbm)
library(rsample)
library(xgboost)
library(caret)
library(h2o)
library(pdp)
library(vip)
library(RColorBrewer)
library(tidyr)

#### IMPORT AND SUBSET DATA ####
### Raw data ###
setwd("~/Google Drive/Cd meta-analysis - Jordon and Andrew")
Cd.dat.all <- read.csv("Data and Stats/Extracted and finalized data.csv")

### Retrieve and add climatic data ###
coord.data <- Cd.dat.all %>% drop_na(Sample.ID, GPS.Lat, GPS.Long)

R <- raster::getData("worldclim", var="bio", res=2.5, download=TRUE)
R <- R[[c(1,12)]]
names(R) <- c("MAT", "MAP")

names <- as.data.frame(coord.data[,1])
coords <- as.data.frame(coord.data[,c(13,12)]) #Long,Lat
coords[,1] <- as.numeric(as.character(coords[,1]))

points <- SpatialPoints(coords, proj4string = R@crs)
values <- raster::extract(R, coords)

clim.dat <- cbind.data.frame(names, coordinates(points), values)
clim.dat$MAT <- (clim.dat$MAT/10)
plot(R[[1]])
plot(points, add=T) # check locations
names(clim.dat)[1] <- "Sample.ID"
Cd.data.all <- merge(x=Cd.dat.all, y=clim.dat, by="Sample.ID", all.x=TRUE)
Cd.data.all <- Cd.data.all[,c(-53,-54)]

Cd.data.all$Country <- as.factor(Cd.data.all$Country)
colnames(Cd.data.all)[which(names(Cd.data.all) == "Country")] <- "region"
Cd.data.all$Study.ID <- as.factor(Cd.data.all$Study.ID)
Cd.data.all$SampleDepth.max <- as.numeric(Cd.data.all$SampleDepth.max)
Cd.data.all$SampleDepth.min <- as.numeric(Cd.data.all$SampleDepth.min)

set.seed(12345)
options(scipen=999)

Cd.data.all %>% dplyr::group_by(region) %>% summarize(mean.pH=mean(Soil.pH, na.rm=TRUE))

#### Figure 1: study data types ####
venn.colors <- brewer.pal(4, "Set3")
d1 <- c("Soil Cd (available)"=7,"Soil Cd (available)&Soil Cd (total)"=4, "Soil Cd (available)&Soil Cd (total)&Bean Cd"=1, "Soil Cd (total)&Bean Cd"=2, "Leaf Cd&Bean Cd"=1, "Soil Cd (available)&Leaf Cd&Bean Cd"=4, "Soil Cd (total)&Leaf Cd&Bean Cd"=2, "Soil Cd (available)&Soil Cd (total)&Leaf Cd&Bean Cd"=2, "Soil Cd (total)"=8)
plot(euler(d1), quantities=TRUE, fills=venn.colors, alpha=0.6, edges=NULL)

plot(venn(d1), fills=venn.colors, alpha=0.6, edges=NULL)
ggsave("Figure 1.pdf", width=5.1, height=4)

#### Figure 2: map and distribution of climate variables ####
### Distribution of samples: country ###
world_map <- subset(map_data("world"), region!="Antarctica")
Cd.by.country <- Cd.data.all %>% dplyr::group_by(region) %>% dplyr::summarize(n.obs=n())
Cd.map.data <- left_join(Cd.by.country, world_map, by="region")
world.map <- ggplot(data=Cd.map.data) + geom_map(data=world_map, map = world_map, aes(map_id=region), fill="white", color="black") + expand_limits(x=world_map$long, y=world_map$lat) + geom_map(map=world_map, aes(map_id=region, fill=n.obs), color="black") + theme_bw() + scale_fill_gradientn(colors=c("#FF0000", "#FBAE00", "#00A58B"), trans="log", breaks=c(10,50,400), labels=c(10,50,400)) + guides(fill=guide_colorbar(direction="horizontal", barwidth=20, barheight=0.6)) + theme(legend.position = "top", legend.title=element_blank()) + ylab("Latitude (degrees)") + xlab("Latitude (degrees)")
world.map
SA.map <- world.map + lims(x=c(-100,-25), y=c(-60,20))
SA.map
Af.map <- world.map + lims(x=c(-25,50), y=c(-30,30)) + theme(legend.position = "none")
Af.map
SEAs.map <- world.map + lims(x=c(90,130), y=c(-30,30)) + theme(legend.position = "none")
SEAs.map

clor.map.all <- ggarrange(SA.map, Af.map, SEAs.map, ncol=3, labels=c("a", "b", "c"), common.legend=TRUE, legend = "top")
clor.map.all


### Climatic data distribution ###
MAT.MAP.plot <- ggplot(data=Cd.data.all, aes(x=MAP, y=MAT)) + geom_point(size=1.6, alpha=0.6) + theme_bw() + labs(x="Mean Annual Precipitation (mm)", y="Mean Annual Temperature (°C)") + theme(legend.position = "bottom",legend.title = element_blank(), legend.text=element_text(size=10), )
MAT.MAP.plot
MAP.MAT.xdens <- axis_canvas(MAT.MAP.plot, axis="x") + geom_density(data=Cd.data.all, aes(x=MAP), alpha=0.6, color=NA, size=0.2, fill="#DC3220")
MAP.MAT.xdens
MAP.MAT.ydens <- axis_canvas(MAT.MAP.plot, axis="y", coord_flip=TRUE) + geom_density(data=Cd.data.all, aes(x=MAT), alpha=0.6, color=NA, size=0.2, fill="#005AB5") + coord_flip()
MAP.MAT.ydens
MAT.MAP.plot1 <- insert_xaxis_grob(MAT.MAP.plot, MAP.MAT.xdens, grid::unit(0.2, "null"), position="top")
MAT.MAP.plot2 <- insert_yaxis_grob(MAT.MAP.plot1, MAP.MAT.ydens, grid::unit(0.2, "null"), position = "right")
MAT.MAP.plot.final <- ggdraw(MAT.MAP.plot2)
MAT.MAP.plot.final

map.n.climate <- ggarrange(clor.map.all, MAT.MAP.plot.final, ncol=1, nrow=2, labels=c("", "d"), heights=c(1,1))
map.n.climate
ggsave("Figure 2.pdf", width=6, height=6)

#### Figure 3: soil depth and Cd concentration ####
# Create fixed country colors #
country_pal <- c("#FF0000", "#00A58B", "#FBAE00", "#FF7F00", "#35BDDB", "#766E66", "#EBCCAE", "#0E6B9A", "#D59D4F", "#ACDEDE", "#000000")
names(country_pal) <- c("Bolivia","Colombia", "Costa Rica", "Ecuador", "Ghana", "Honduras", "Indonesia","Malaysia","Nigeria", "Peru", "Trinidad and Tobago")
print(country_pal)

### Available Cd ###
avail.data <- Cd.data.all %>% drop_na(Soil.AvailCd1) %>% droplevels()

avail.depth.count <- ggplot(data=avail.data, aes(x=SampleDepth.max)) + geom_histogram(binwidth = 5, fill="#AFABAB") + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + ylab("Count") + theme_bw() + coord_flip(xlim=c(140,0)) + theme(legend.position = "none")
avail.depth.count

avail.Cd.depth <- ggplot(data=avail.data, aes(x=SampleDepth.max, y=Soil.AvailCd1)) + geom_jitter(aes(color=region), alpha=0.8, size=0.75) + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + labs(y=expression("Available Soil Cd (mg"~kg^-1~"soil)")) + geom_smooth(method="loess", se=TRUE, level=0.95, color="#424242") + theme_bw() + coord_flip(ylim=c(0,5), xlim=c(140,0)) + scale_color_manual(values=country_pal) + theme(legend.position = "none")
avail.Cd.depth

avail.depth <- ggarrange(avail.depth.count, avail.Cd.depth, ncol=2, labels=c("a", "b"), align="h", heights=c(2,2)) # will be combined with totals
avail.depth

### Total Cd ###
total.data <- Cd.data.all %>% drop_na(Soil.TotalCd) %>% droplevels()
total.data.NoArguello <- subset(total.data, LeadAuthor.primary !="David Arguello")

total.depth.count <- ggplot(data=total.data, aes(x=SampleDepth.max)) + geom_histogram(binwidth = 5, fill="#AFABAB") + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + ylab("Count") + theme_bw() + coord_flip(xlim=c(140,0)) + theme(legend.position = "none")
total.depth.count

total.depth.count.noArg <- ggplot(data=total.data.NoArguello, aes(x=SampleDepth.max)) + geom_histogram(binwidth = 5, fill="#AFABAB") + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + ylab("Count") + theme_bw() + coord_flip(xlim=c(140,0))
total.depth.count.noArg

total.Cd.depth <- ggplot(data=total.data, aes(x=SampleDepth.max, y=Soil.TotalCd)) + geom_jitter(aes(color=region), alpha=0.8, size=0.75) + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + labs(y=expression("Total Soil Cd (mg"~kg^-1~"soil)")) + geom_smooth(method="loess", level=0.95, se=TRUE, color="#424242") + scale_color_manual(values=country_pal) + theme_bw() + coord_flip(ylim=c(0,20), xlim=c(140,0)) + theme(legend.position="none")
total.Cd.depth

total.Cd.depth.noArg <- ggplot(data=total.data.NoArguello, aes(x=SampleDepth.max, y=Soil.TotalCd)) + geom_jitter(aes(color=region), alpha=0.8, size=0.75) + scale_x_reverse(breaks=seq(150,0, -30)) + xlab("Sample depth (cm)") + labs(y=expression("Total Soil Cd (mg Cd"~kg^-1~"soil)")) + geom_smooth(method="loess", se=TRUE, level=0.95, color="#424242") + scale_color_manual(values=country_pal) + theme_bw() + coord_flip(ylim=c(0,20), xlim=c(140,0)) + guides(color=guide_legend(override.aes=list(size=1))) + theme(legend.title = element_blank())
total.Cd.depth.noArg

total.depth <- ggarrange(avail.depth.count, avail.Cd.depth, total.depth.count, total.Cd.depth, total.depth.count.noArg, total.Cd.depth.noArg, ncol=2, nrow=3, labels=c("a","b", "c", "d", "e", "f"), align="hv", widths=c(1.5,2,1.5,2,1.5,2), heights=c(2,2,2,2,2,2), common.legend = TRUE, legend="top")
total.depth
ggsave("Figure 3 (prelim).pdf", width=6, height=8.5)

#### BRT MODELS: whole bean Cd ####
set.seed(12345)
### Wrangle data to get depth-weighted soil Cd ###
gbm.dat <- as.data.frame(Cd.data.all %>% dplyr::select(Study.ID, SiteName, Block.Rep, SampleDepth.min, SampleDepth.max, Soil.TotalCd, Soil.AvailCd1, WholeBeanCd, LeafCd, GPS.Lat.x, Clay.pct, Soil.pH, CEC, SOC.pct, MAT, MAP) %>% dplyr::group_by(Study.ID, SiteName) %>% mutate(layer_depth = (SampleDepth.max - SampleDepth.min), site_weight = layer_depth/max(SampleDepth.max), TotalCd.wtd = site_weight*Soil.TotalCd, AvailCd.wtd = site_weight*Soil.AvailCd1))

gbm.data <- as.data.frame(gbm.dat %>% dplyr::group_by(Study.ID, SiteName, GPS.Lat.x, Clay.pct, Soil.pH, CEC, SOC.pct, MAT, MAP, WholeBeanCd, LeafCd) %>% summarize(TotalCd.wtd = sum(TotalCd.wtd), AvailCd.wtd =sum(AvailCd.wtd)) %>% mutate(BCFt = WholeBeanCd/TotalCd.wtd, BCFa = WholeBeanCd/AvailCd.wtd) %>% drop_na(WholeBeanCd))

### Do duh model ###
gbm.split <- initial_split(gbm.data, prop=0.8)
gbm.train <- training(gbm.split)
gbm.test <- testing(gbm.split)

## Tuning Model ##
hyper_grid <- expand.grid(shrinkage = c(.001, .05, .075, .1, .12, .14, .16, .18, .2, .22, .24, .26, .28, .3, .325, .35, .4), interaction.depth = c(1,2,3,4,5,6,7,8,9), n.minobsinnode = c(5,6,7,8,9,10,11,12,13,14,15), bag.fraction = c(.65,.7,.75,.8,.85,.9), optimal_trees = 0, min_RMSE = 0)

nrow(hyper_grid)

random_index <- sample(1:nrow(gbm.train), nrow(gbm.train))
random_gbm_train <- gbm.train[random_index, ]

pb <- txtProgressBar(min = 0, max = 10098, style = 3); for(i in 1:nrow(hyper_grid)) {
  
  gbm.tune <- gbm(formula = WholeBeanCd ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd + AvailCd.wtd + TotalCd.wtd, distribution = "gaussian", data = random_gbm_train, n.trees = 2000, interaction.depth = hyper_grid$interaction.depth[i], shrinkage = hyper_grid$shrinkage[i], n.minobsinnode = hyper_grid$n.minobsinnode[i], bag.fraction = hyper_grid$bag.fraction[i], train.fraction = .75, n.cores = NULL, verbose = FALSE)
  
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
  
  setTxtProgressBar(pb, i)
}; close(pb)

hyper_grid %>% dplyr::arrange(min_RMSE) %>% head(10)

## Finalized model ##
set.seed(12345)
gbm.fit.final <- gbm(formula=WholeBeanCd ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + TotalCd.wtd + AvailCd.wtd + LeafCd, distribution="gaussian", data=gbm.train, n.trees=11, interaction.depth = 9, shrinkage=0.35, n.minobsinnode = 12, bag.fraction = 0.65, train.fraction=1)

print(gbm.fit.final)
summary(gbm.fit.final, cBars=10, method=permutation.test.gbm, las=2, normalize=TRUE)

summary(gbm.data$Soil.pH)

#### Figure 4: bean Cd drivers and PDP plots ####
VIP <- c(50.5, 24.7, 11.9, 7.6, 3.9, 0.8, 0.6, 0.0, 0.0)
Variable <- c("Total Soil Cd", "Soil pH", "Leaf Cd", "SOC", "CEC", "MAP", "Clay", "Available Soil Cd", "MAT")
sig <- c("Y", "Y", "Y", "N", "N", "N", "N", "N", "N")
VIP.data <- data.frame(Variable, VIP, sig, stringsAsFactors = FALSE)
VIP.data$Variable <- as.factor(VIP.data$Variable)
VIP.data$sig <- as.factor(VIP.data$sig)
str(VIP.data)
VIP.dot <- ggdotchart(data=VIP.data, x="Variable", y="VIP", color="sig", palette=c("#999999","#DC3220"), sorting="descending", rotate=TRUE, add="segments", add.params=list(color="lightgrey", size=2), group="sig", dot.size=7, label=formatC(round(VIP.data$VIP,1), format="f", digits=1), font.label=list(color = "black", size = 7, vjust = 0.5), ggtheme=theme_bw()) + theme(axis.title.y=element_blank(), legend.position="none")
VIP.dot

# plot partial dependency plots #
Total.Cd.gbm.plot <- gbm.fit.final %>% partial(pred.var = "TotalCd.wtd", n.trees = gbm.fit.final$n.trees, grid.resolution = 150) %>% autoplot(train = gbm.train, smooth=TRUE, smooth.method="loess", smooth.span=0.27, lty="dashed") + theme_bw() + labs(y="Bean Cd\n(mg Cd / kg dry matter)", x="Depth-weighted total soil Cd\n(mg Cd / kg soil)") + xlim(0,4)
Total.Cd.gbm.plot

pH.gbm.plot <- gbm.fit.final %>% partial(pred.var = "Soil.pH", n.trees = gbm.fit.final$n.trees, grid.resolution = 100) %>% autoplot(train = gbm.train, smooth=TRUE, smooth.method="loess", smooth.span=0.37, lty="dashed") + theme_bw() + ylab("Bean Cd\n(mg Cd / kg dry matter)") + xlab("Soil pH") + xlim(4,7)
pH.gbm.plot

Leaf.Cd.gbm.plot <- gbm.fit.final %>% partial(pred.var = "LeafCd", n.trees = gbm.fit.final$n.trees, grid.resolution = 150) %>% autoplot(train = gbm.train, smooth=TRUE, smooth.method="loess", smooth.span=0.27, lty="dashed") + theme_bw() + labs(y="Bean Cd\n(mg Cd / kg dry matter)", x="Leaf Cd\n(mg Cd / kg dry matter)") + xlim(0,8)
Leaf.Cd.gbm.plot

gbm.plots <- ggarrange(VIP.dot, Total.Cd.gbm.plot, pH.gbm.plot, Leaf.Cd.gbm.plot, ncol=2, nrow=2, labels=c("a","b", "c", "d"), align="h")
gbm.plots
ggsave("Figure 4.pdf", width=8, height=6)

pred.bean <- predict(gbm.fit.final, n.trees = gbm.fit.final$n.trees, gbm.test)
caret::RMSE(pred.bean, gbm.test$WholeBeanCd)

#### BRT MODELS: bioconcentration factor ####
### Make-uh duh model ###
set.seed(12345)
BCF.tot.data <- as.data.frame(gbm.data %>% filter(TotalCd.wtd > 0) %>% drop_na(BCFt))

BCFt.split <- initial_split(BCF.tot.data, prop=0.8)
BCFt.train <- training(BCFt.split)
BCFt.test <- testing(BCFt.split)

## TUNING ##
BCF_hyper_grid <- expand.grid(shrinkage = c(.001, .05, .075, .1, .12, .14, .16, .18, .2, .22, .24, .26, .28, .3, .325, .35, .4), interaction.depth = c(1,2,3,4,5,6,7,8,9), n.minobsinnode = c(5,6,7,8,9,10,11,12,13,14,15), bag.fraction = c(.65,.7,.75,.8,.85,.9), optimal_trees = 0, min_RMSE = 0)

nrow(BCF_hyper_grid)

# BCF only #
random_index_t <- sample(1:nrow(BCFt.train), nrow(BCFt.train))
random_BCFt_train <- BCFt.train[random_index_t, ]

pb_BCF <- txtProgressBar(min = 0, max = 10098, style = 3); for(i in 1:nrow(BCF_hyper_grid)) {
  
  BCF.tune <- gbm(formula = BCFt ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd + TotalCd.wtd, distribution = "gaussian", data = random_BCFt_train, n.trees = 2000, interaction.depth = BCF_hyper_grid$interaction.depth[i], shrinkage = BCF_hyper_grid$shrinkage[i], n.minobsinnode = BCF_hyper_grid$n.minobsinnode[i], bag.fraction = BCF_hyper_grid$bag.fraction[i], train.fraction = .75, n.cores = NULL, verbose = FALSE)
  
  BCF_hyper_grid$optimal_trees[i] <- which.min(BCF.tune$valid.error)
  BCF_hyper_grid$min_RMSE[i] <- sqrt(min(BCF.tune$valid.error))
  
  setTxtProgressBar(pb_BCF, i)
}; close(pb_BCF)

BCF_hyper_grid %>% dplyr::arrange(min_RMSE) %>% head(10)

# Fit the final model #
BCFt.fit.final <- gbm(formula= BCFt ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd + TotalCd.wtd, distribution="gaussian", data=BCFt.train, n.trees=17, interaction.depth = 4, shrinkage=0.40, n.minobsinnode = 10, bag.fraction = 0.75, train.fraction=1)

print(BCFt.fit.final)
summary(BCFt.fit.final, cBars=10, method=permutation.test.gbm, las=2, normalize=TRUE)

pred.BCF <- predict(BCFt.fit.final, n.trees = BCFt.fit.final$n.trees, BCFt.test)
caret::RMSE(pred.BCF, BCFt.test$WholeBeanCd)

## Plots of distribution, VIP, and PDPs ##
# Distribution of values #
summary(BCF.tot.data$BCFt)
sum(BCF.tot.data$BCFt < 1) / length(BCF.tot.data$BCFt)
sum(BCF.tot.data$BCFt > 1 & BCF.tot.data$BCFt < 5) / length(BCF.tot.data$BCFt)
sum(BCF.tot.data$BCFt > 5) / length(BCF.tot.data$BCFt)
quantile(BCF.tot.data$BCFt, c(0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95))

#### Figure 5: BCF plots ####
# density plot #
BCFt.densityplot <- ggplot(data=BCF.tot.data, aes(x=BCFt)) + geom_density(fill="#005AB5", color=NA) + theme_bw()
BCFt.d <- ggplot_build(BCFt.densityplot)$data[[1]]
BCFt.plot <- BCFt.densityplot + geom_area(data=subset(BCFt.d, x > 4.9), aes(x=x, y=y), fill="#DC3220") + geom_area(data=subset(BCFt.d, x > 1 & x < 5), aes(x=x, y=y), fill="#FFBC42") + labs(y="Density", x="Bioconcentration Factor (Bean Cd : Total Soil Cd)")
BCFt.plot

# Soil pH #
BCFt.pH.plot <- BCFt.fit.final %>% partial(pred.var = "Soil.pH", n.trees = BCFt.fit.final$n.trees, grid.resolution = 150) %>% autoplot(train = BCFt.train, smooth=TRUE, smooth.method="loess", smooth.span=0.32, lty="dashed") + theme_bw() + labs(y="Bioconcentration Factor", x="Soil pH") + xlim(4,7)
BCFt.pH.plot

# SOC #
BCFt.SOC.plot <- BCFt.fit.final %>% partial(pred.var = "SOC.pct", n.trees = BCFt.fit.final$n.trees, grid.resolution = 200) %>% autoplot(train = BCFt.train, smooth=TRUE, smooth.method="loess", smooth.span=0.32, lty="dashed") + theme_bw() + labs(y="Bioconcentration Factor", x="SOC (%)") + xlim(0,10)
BCFt.SOC.plot

BCFt.bottom.plots <- ggarrange(BCFt.pH.plot, BCFt.SOC.plot, ncol=2, nrow=1, labels=c("b", "c"), align="v")
BCFt.bottom.plots

BCFt.plots.all <- ggarrange(BCFt.plot, BCFt.bottom.plots, ncol=1, nrow=2, labels=c("a", ""), heights=c(0.8,1))
BCFt.plots.all
ggsave("Figure 5 (prelim).pdf", width=8, height=6)

#### SUPPLEMENTARY INFORMATION THINGS ####
### Figure S2: Correlation matrix ###
Cd.cor.data <- Cd.data.all[,c(42,26,24,53,54,17,21:23)]
names(Cd.cor.data)[names(Cd.cor.data) == "SOC.pct"] <- "SOC"
names(Cd.cor.data)[names(Cd.cor.data) == "Soil.pH"] <- "Soil pH"
names(Cd.cor.data)[names(Cd.cor.data) == "Clay.pct"] <- "Clay"
names(Cd.cor.data)[names(Cd.cor.data) == "Soil.TotalCd"] <- "Total Soil Cd"
names(Cd.cor.data)[names(Cd.cor.data) == "Soil.AvailCd1"] <- "Avail. Soil Cd"
names(Cd.cor.data)[names(Cd.cor.data) == "LeafCd"] <- "Leaf Cd"

cormat <- round(cor(Cd.cor.data, use="pairwise.complete.obs"),2)
cormat
melt_cormat <- melt(cormat)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
melt_cormat <- melt(upper_tri, na.rm = TRUE)

cormat.plot <- ggplot(data = melt_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") + scale_fill_gradient2(low = "#005AB5", high = "#DC3220", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation\nCoefficient") + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(), legend.justification = c(1, 0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
cormat.plot
ggsave("Figure S2.pdf", width=5, height=4.7)

### Robustness checks ###
## Bean Cd ##
# Model without leaf Cd #
gbm.fit.NoLeaf <- gbm(formula=WholeBeanCd ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + TotalCd.wtd + AvailCd.wtd, distribution="gaussian", data=gbm.train, n.trees=23, interaction.depth = 5, shrinkage=0.35, n.minobsinnode = 5, bag.fraction = 0.75, train.fraction=1)
print(gbm.fit.NoLeaf)
summary(gbm.fit.NoLeaf, cBars=10, method=permutation.test.gbm, las=2, normalize=TRUE)

pred.nl <- predict(gbm.fit.NoLeaf, n.trees = gbm.fit.NoLeaf$n.trees, gbm.test)
caret::RMSE(pred.nl, gbm.test$WholeBeanCd)

# Model without total soil Cd #
gbm.fit.NoTot <- gbm(formula=WholeBeanCd ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd + AvailCd.wtd, distribution="gaussian", data=gbm.train, n.trees=3, interaction.depth = 5, shrinkage=0.40, n.minobsinnode = 7, bag.fraction = 0.70, train.fraction=1)
print(gbm.fit.NoTot)
summary(gbm.fit.NoTot, cBars=10, method=permutation.test.gbm, las=2, normalize=TRUE)

pred.nt <- predict(gbm.fit.NoTot, n.trees = gbm.fit.NoTot$n.trees, gbm.test)
caret::RMSE(pred.nt, gbm.test$WholeBeanCd)

# Individual conditional expectation (ICE) plots #
ice.TotCd <- gbm.fit.final %>%
  partial(
    pred.var = "TotalCd.wtd", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 150,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = gbm.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("Total Soil Cd (mg Cd / kg soil)") + xlim(0,4) + theme_bw()
ice.TotCd

ice.pH <- gbm.fit.final %>%
  partial(
    pred.var = "Soil.pH", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 100,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = gbm.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("Soil pH")+ xlim(4,7) + theme_bw()
ice.pH

ice.LeafCd <- gbm.fit.final %>%
  partial(
    pred.var = "LeafCd", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 150,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = gbm.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("Leaf Cd (mg Cd / kg mass)") + xlim(0,8) + theme_bw()
ice.LeafCd

ICE.plots <- ggarrange(ice.TotCd, ice.pH, ice.LeafCd, ncol=1, nrow=3, labels=c("a","b", "c"), align="v")
ICE.plots
ggsave("Figure S3.pdf", width=5, height=8)

SOC.gbm.plot <- gbm.fit.final %>% partial(pred.var = "SOC.pct", n.trees = gbm.fit.final$n.trees, grid.resolution = 150) %>% autoplot(train = gbm.train, smooth=TRUE, smooth.method="loess", smooth.span=0.2, lty="dashed") + theme_bw() + labs(y="Bean Cd\n(mg Cd / kg dry matter)", x="SOC (%)")
SOC.gbm.plot

ice.SOC <- gbm.fit.final %>%
  partial(
    pred.var = "SOC.pct", 
    n.trees = gbm.fit.final$n.trees, 
    grid.resolution = 150,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = gbm.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("SOC (%)") + theme_bw()
ice.SOC
SOC.plots <- ggarrange(SOC.gbm.plot, ice.SOC, ncol=1, nrow=2, labels=c("a","b"), align="v")
SOC.plots
ggsave("Figure S4.pdf", width=4, height=6)

## Bioconcentration factor ##
# Model without total soil Cd #
BCF_hyper_grid <- expand.grid(shrinkage = c(.001, .05, .075, .1, .12, .14, .16, .18, .2, .22, .24, .26, .28, .3, .325, .35, .4), interaction.depth = c(1,2,3,4,5,6,7,8,9), n.minobsinnode = c(5,6,7,8,9,10,11,12,13,14,15), bag.fraction = c(.65,.7,.75,.8,.85,.9), optimal_trees = 0, min_RMSE = 0)

nrow(BCF_hyper_grid)

# BCF only #
set.seed(12345)
random_index_t <- sample(1:nrow(BCFt.train), nrow(BCFt.train))
random_BCFt_train <- BCFt.train[random_index_t, ]

pb_BCF_noT <- txtProgressBar(min = 0, max = 10098, style = 3); for(i in 1:nrow(BCF_hyper_grid)) {
  
  BCF.tune.noT <- gbm(formula = BCFt ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd, distribution = "gaussian", data = random_BCFt_train, n.trees = 2000, interaction.depth = BCF_hyper_grid$interaction.depth[i], shrinkage = BCF_hyper_grid$shrinkage[i], n.minobsinnode = BCF_hyper_grid$n.minobsinnode[i], bag.fraction = BCF_hyper_grid$bag.fraction[i], train.fraction = .75, n.cores = NULL, verbose = FALSE)
  
  BCF_hyper_grid$optimal_trees[i] <- which.min(BCF.tune.noT$valid.error)
  BCF_hyper_grid$min_RMSE[i] <- sqrt(min(BCF.tune.noT$valid.error))
  
  setTxtProgressBar(pb_BCF_noT, i)
}; close(pb_BCF_noT)

BCF_hyper_grid %>% dplyr::arrange(min_RMSE) %>% head(10)

# Fit the final model #
BCFt.fit.final.noT <- gbm(formula= BCFt ~ Clay.pct + Soil.pH + CEC + SOC.pct + MAT + MAP + LeafCd, distribution="gaussian", data=BCFt.train, n.trees=5, interaction.depth = 9, shrinkage=0.40, n.minobsinnode = 11, bag.fraction = 0.80, train.fraction=1)

print(BCFt.fit.final.noT)
summary(BCFt.fit.final.noT, cBars=10, method=permutation.test.gbm, las=2, normalize=TRUE)

pred.BCF.noT <- predict(BCFt.fit.final.noT, n.trees = BCFt.fit.final.noT$n.trees, BCFt.test)
caret::RMSE(pred.BCF.noT, BCFt.test$WholeBeanCd)

## Plots of distribution, VIP, and PDPs ##
# Final model: Individual conditional expectation (ICE) plots #
BCFt.ice.pH <- BCFt.fit.final %>%
  partial(
    pred.var = "Soil.pH", 
    n.trees = BCFt.fit.final$n.trees, 
    grid.resolution = 150,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = BCFt.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("Soil pH") + theme_bw()
BCFt.ice.pH

BCFt.ice.SOC <- BCFt.fit.final %>%
  partial(
    pred.var = "SOC.pct", 
    n.trees = BCFt.fit.final$n.trees, 
    grid.resolution = 150,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = BCFt.train, alpha = .1, center=TRUE) + ylab("Change in Bean Cd \n(mg Cd / kg mass)") + xlab("SOC (%)") + theme_bw()
BCFt.ice.SOC

BCFt.ice.plots <- ggarrange(BCFt.ice.pH, BCFt.ice.SOC, ncol=1, nrow=2, labels=c("a", "b"), align="v")
BCFt.ice.plots
ggsave("Figure S5.pdf", height=6, width=4)