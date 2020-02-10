# CART model to compare conditions where model simulations under/over predict vector data ------------
rm(list=ls()) #remove previous variable assignments

# load libraries 
library(rpart)
library(rattle)

# mosquitoes
rpart.tree.adults <- rpart(Adult_correspondence_magnitude ~ T_min + T_mean + T_max + T_var + H_min + H_mean + H_max + H_var + R_min + R_mean + R_max + R_var + Site + Country + Urbanization + CoastalInland, data=zscore_aedes)
fancyRpartPlot(rpart.tree.adults)
prune.rpart.tree.adults <- prune(rpart.tree.adults, cp=0.04) # pruning the tree
fancyRpartPlot(prune.rpart.tree.adults)

# arboviruses
rpart.tree.arbo <- rpart(Arboviruses_correspondence_magnitude ~ T_min + T_mean + T_max + T_var + H_min + H_mean + H_max + H_var + R_min + R_mean + R_max + R_var + Country + Urbanization + CoastalInland, data=zscore_arboviruses)
fancyRpartPlot(rpart.tree.arbo)
prune.rpart.tree.arbo <- prune(rpart.tree.arbo, cp=0.03) # pruning the tree
fancyRpartPlot(prune.rpart.tree.arbo)

