# Visually compare model output given different initial conditions ------------------------
rm(list=ls()) #remove previous variable assignments

# load libraries 
library(ggplot2)

# load data 
traits <- read.csv("output/SEI-SEIR_simulations_with_trait_variation.csv", head=T)
load("output/CCF.RData")

# adjust dates and ic classes
traits$Date <- as.Date(traits$Date, "%Y-%m-%d")
traits$simulation_number <- as.factor(traits$simulation_number)

traits2 <- subset(traits, time > 90)
traits2 <- subset(traits2, Date > "2014-01-01") 

# plot trait sensitvitity analysis for cases
cases <- subset(df_max_ccf, Y=="cases")
cases <- merge(cases, traits2, by=c("Site", "Rain_function"))

ggplot() +
  geom_line(data = cases, aes(x = Date, y = I, color = as.factor(simulation_number), group = as.factor(simulation_number))) +
  facet_wrap(~Site, ncol=2, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Predicted arboviral cases")

# plot trait sensitvitity analysis for vectors
vectors <- subset(df_max_ccf, Y=="vectors")
vectors <- merge(vectors, traits2, by=c("Site", "Rain_function"))
vectors$Mtot <- vectors$M1 + vectors$M2 + vectors$M3

ggplot() +
  geom_line(data = vectors, aes(x = Date, y = Mtot, color = as.factor(simulation_number), group = as.factor(simulation_number))) +
  facet_wrap(~Site, ncol=2, scales = "free") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Predicted mosquito abundance")
