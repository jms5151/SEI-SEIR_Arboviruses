# Temperature, humidity, and rainfall dependent SEI-SEIR model simulations -----------------------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load data 
source("SEI-SEIR_model_TRH.R")
source("SEI-SEIR_simulation_setup.R")

# run simulations
traitDF <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "Site")
traitFileName <- "SEI-SEIR_simulations_TRH_final_model.csv"
write.csv(traitDF, traitFileName, row.names = F)

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
  out2$Site <- sites[l]
  traitDF <- rbind(traitDF, out2)
  write.csv(traitDF, traitFileName, row.names = F)
  cat("finished running ode for", sites[l], "\n")
}
