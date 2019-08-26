# SEI-SEIR model simulations with trait variability --------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load model
source("Codes/SEI-SEIR_model_with_trait_variation.R")

# load data 
source("Codes/SEI-SEIR_simulation_setup.R")

# run simulations
traitDF <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "simulation_number", "Site")
traitFileName <- "Concatenated_Data/model_simulations/SEI-SEIR_simulations_with_trait_variation.csv"
write.csv(traitDF, traitFileName, row.names = F)

for (j in 1:nrow(trait_posterior)){
  for (k in 2:ncol(trait_posterior)){
    x <- colnames(trait_posterior)[k]
    y <- trait_posterior[j,k]
    assign(x,y)
  }
  for (l in 1:length(sites)){
    climateData2 <- subset(climateData, Site == sites[l])
    climateData2 <- climateData2[order(climateData2$Date),]
    climateData2 <- climateData2[complete.cases(climateData2),]
    temp <- climateData2$Temperature
    rain <- climateData2$Two_week_rainfall
    Rmax <- 123
    if (unique(climateData2$country) == "Ecuador"){
      K_trh <- K_trh_quadratic
    } else {
      K_trh <- K_trh_briere
    }
    hum <- climateData2$SVPD
    Date <- climateData2$Date
    H0 <- population[l]
    city <- sites[l]
    BR <- BRs[l]
    DR <- DRs[l]
    times <- seq(1,length(Date), by=1)
    M0 <- K_trh(temp[1], hum[1], mean(rain), Rmax, H0)
    parameters <- c(EFD, pEA, MDR, K_trh, a, pMI, mu_th, PDR, b, timestep=timestep)
    state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*H0, E = startIC$e*H0, I = startIC$i*H0, R = startIC$r*H0)
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$simulation_number <- j
    out2$Site <- sites[l]
    write.table(out2, file=traitFileName, row.names = F, sep = ",", col.names = !file.exists(traitFileName), append = T)
    cat("finished running ode for", sites[l], "simulation #", j, "\n")
  }
}
