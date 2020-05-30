# SEI-SEIR model simulations with trait variability --------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load model
source("models/SEI-SEIR_model_with_trait_variation.R")

# load data 
source("models/SEI-SEIR_simulation_setup.R")
load("output/CCF.RData")
df_max_ccf <- df_max_ccf[!duplicated(df_max_ccf[c("Site", "Rain_function")]),]

# set population
pop <- data.frame("Site" = sites, "population" = population, "BR" = BRs, "DR" = DRs)
df_max_ccf <- merge(df_max_ccf, pop, by="Site")

# run simulations
traitDF <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "simulation_number", "Site", "Rain_function")
traitFileName <- "output/SEI-SEIR_simulations_with_trait_variation.csv"
write.csv(traitDF, traitFileName, row.names = F)

rfunctions_names <- c("Briere", "Quadratic", "Inverse")
rfunctions <- list(K_thr_briere, K_thr_quadratic, K_thr_inverse)

for (j in 1:nrow(trait_posterior)){
  for (k in 2:ncol(trait_posterior)){
    x <- colnames(trait_posterior)[k]
    y <- trait_posterior[j,k]
    assign(x,y)
  }
  for (l in 1:nrow(df_max_ccf)){
    climateData2 <- subset(climateData, Site == df_max_ccf$Site[l])
    climateData2 <- climateData2[order(climateData2$Date),]
    climateData2 <- climateData2[complete.cases(climateData2),]
    temp <- climateData2$Temperature
    rain <- climateData2$Two_week_rainfall
    Rmax <- 123
    rf <- which(rfunctions_names==df_max_ccf$Rain_function[l])
    K_thr <- rfunctions[[rf]]
    hum <- climateData2$SVPD
    Date <- climateData2$Date
    N <- df_max_ccf$population[l]
    city <- df_max_ccf$Site[l]
    BR <- df_max_ccf$BR[l]
    DR <- df_max_ccf$DR[l]
    times <- seq(1,length(Date), by=1)
    M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
    state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$simulation_number <- j
    out2$Site <- df_max_ccf$Site[l]
    out2$Rain_function <- df_max_ccf$Rain_function[l]
    write.table(out2, file=traitFileName, row.names = F, sep = ",", col.names = !file.exists(traitFileName), append = T)
    cat("finished running ode for", df_max_ccf$Site[l], "simulation #", j, "\n")
  }
}
