# data to set up model simulations -------------------------------------------------------------
# load packages
library(deSolve)

# load climate data
load("Concatenated_Data/climate_data/merged_climate_data.RData")

# load and set initial conditions
init.cond <- read.csv("LHS_inputs.csv", head=T)
startIC <- subset(init.cond, IC == "18")

# load traits
trait_posterior <- read.csv("Random_sample_of_posterior_traits.csv", head=T)

# set immigration and emmigration rate
ie <- 0.01

# set up list of sites
sites <- c("Chulaimbo", "Kisumu", "Msambweni", "Ukunda", "Huaquillas", "Machala", "Portovelo", "Zaruma")

# set human population numbers for each site 
population <- c(7304, 547557, 240698, 154048, 57370, 279890, 13670, 25620)

# set birth and death rates
BRs <- c(rep(31.782,4),rep(20.175,4)) # birth rates from https://data.worldbank.org/indicator/SP.DYN.CBRT.IN
DRs <- c(rep(5.284,4),rep(5.121,4)) # death rates from https://data.worldbank.org/indicator/SP.DYN.CBRT.IN

# model timestep
timestep = 1/12