# data to set up model simulations -------------------------------------------------------------
# load packages
library(deSolve)

# load climate data
load("data/climate_data.RData")

# load and set initial conditions
load("data/LHS_inputs.RData")
startIC <- subset(init.cond, IC == "18")

# load traits
load("data/Random_sample_of_posterior_traits.RData")

# set immigration and emmigration rate
ie <- 0.01

# set up list of sites
sites <- c("Chulaimbo", "Kisumu", "Msambweni", "Ukunda", "Huaquillas", "Machala", "Portovelo", "Zaruma")

# set human population numbers for each site 
population <- c(7304, 547557, 240698, 154048, 57366, 279887, 13673, 25615)

# set birth and death rates
BRs <- c(rep(31.782,4),rep(20.175,4)) # birth rates from https://data.worldbank.org/indicator/SP.DYN.CBRT.IN
DRs <- c(rep(5.284,4),rep(5.121,4)) # death rates from https://data.worldbank.org/indicator/SP.DYN.CBRT.IN

# model timestep
timestep = 1/12