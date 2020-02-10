# Run SEI-SEIR model with different rainfall functions ------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load data 
source("models/SEI-SEIR_model_THR.R")
source("models/SEI-SEIR_simulation_setup.R")

# run simulations
traitDF <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "Site", "Rain_function")
traitFileName <- "SEI-SEIR_simulations_THR.csv"
write.csv(traitDF, traitFileName, row.names = F)

rfunctions_names <- c("Briere", "Quadratic", "Inverse")
rfunctions <- list(K_thr_briere, K_thr_quadratic, K_thr_inverse)

for (i in 1:length(sites)){
  climateData2 <- subset(climateData, Site == sites[i])
  temp <- climateData2$Temperature
  rain <- climateData2$Two_week_rainfall
  Rmax <- 123
  hum <- climateData2$SVPD
  Date <- climateData2$Date
  H0 <- population[i]
  city <- sites[i]
  BR <- BRs[i]
  DR <- DRs[i]
  times <- seq(1,length(Date), by=1)
  for (k in 1:length(rfunctions_names)){
    K_thr <- rfunctions[[k]]
    M0 <- K_thr(temp[1], mean(rain), Rmax, H0, timestep)
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
    state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*H0, E = startIC$e*H0, I = startIC$i*H0, R = startIC$r*H0)
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$Site <- sites[i]
    out2$Rain_function <- rfunctions_names[k]
    traitDF <- rbind(traitDF, out2)
    write.csv(traitDF, traitFileName, row.names = F)
    cat("finished running ode for", sites[i], rfunctions_names[k], "\n")
  }
}
