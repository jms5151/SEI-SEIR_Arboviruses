# CART model to compare conditions where model simulations under/over predict vector data ------------
rm(list=ls()) #remove previous variable assignments

# load libraries -------------------------------------------------------------------------------------
library(rpart)
library(rattle)

# load data ------------------------------------------------------------------------------------------
load("zscores_with_climate_data.RData")
zscoreSub$Country <- ifelse(zscoreSub$Site=="Chulaimbo"|zscoreSub$Site=="Kisumu"|zscoreSub$Site=="Msambweni"|zscoreSub$Site=="Ukunda", "Kenya", "Ecuador")
zscoreSub$Urbanization <- ifelse(zscoreSub$Site=="Chulaimbo"|zscoreSub$Site=="Msambweni"|zscoreSub$Site=="Zaruma"|zscoreSub$Site=="Portovelo", "Rural", "Urban")
zscoreSub$CoastalInland <- ifelse(zscoreSub$Site=="Chulaimbo"|zscoreSub$Site=="Kisumu"|zscoreSub$Site=="Zaruma"|zscoreSub$Site=="Portovelo"|zscoreSub$Site=="Huaquillas", "Inland", "Coastal")

# mosquitoes
rpart.tree.adults <- rpart(Adult_correspondence_magnitude ~ T_min + T_mean + T_max + T_var + H_min + H_mean + H_max + H_var + R_min + R_mean + R_max + R_var + Site + Country + Urbanization + CoastalInland, data=zscoreSub)
fancyRpartPlot(rpart.tree.adults)
prune.rpart.tree.adults <- prune(rpart.tree.adults, cp=0.05) # pruning the tree
fancyRpartPlot(prune.rpart.tree.adults)

# dengue
rpart.tree.dengue <- rpart(Dengue_correspondence_magnitude ~ T_min + T_mean + T_max + T_var + H_min + H_mean + H_max + H_var + R_min + R_mean + R_max + R_var + Site + Country + Urbanization + CoastalInland, data=zscoreSub)
fancyRpartPlot(rpart.tree.dengue)
prune.rpart.tree.dengue <- prune(rpart.tree.dengue, cp=0.035) # pruning the tree
fancyRpartPlot(prune.rpart.tree.dengue)

