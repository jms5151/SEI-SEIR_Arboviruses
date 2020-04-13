# SEI-SEIR model simulations with intervention: reduce contact rate ---------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load data 
source("Codes/SEI-SEIR_model_with_intervention_reduceBiteRate.R")
source("Codes/SEI-SEIR_simulation_setup.R")

# create intervention strategy
reductionStrategies <- c(0.00, 0.10, 0.50, 0.90)

# run simulations
traitDF <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "Site", "Intervention_strategy")
traitFileName <- "Concatenated_Data/model_simulations/SEI-SEIR_simulations_with_intervention_reduceBiteRate.csv"
write.csv(traitDF, traitFileName, row.names = F)

for (l in 1:length(sites)){
  climateData2 <- subset(climateData, Site == sites[l] & Date > "2016-10-01" & Date < "2018-03-31")
  climateData2 <- climateData2[order(climateData2$Date),]
  climateData2 <- climateData2[complete.cases(climateData2),]
  temp <- climateData2$Temperature
  rain <- climateData2$Two_week_rainfall
  Rmax <- 123
  if (sites[l] == "Msambweni"){
    K_thr <- K_thr_inverse
  } else if (sites[l] == "Machala"){
    K_thr <- K_thr_quadratic
  } else {
    K_thr <- K_thr_briere
  }
  hum <- climateData2$SVPD
  Date <- climateData2$Date
  N <- population[l]
  city <- sites[l]
  BR <- BRs[l]
  DR <- DRs[l]
  times <- seq(1,length(Date), by=1)
  M0 <- K_thr(temp[1], mean(rain), Rmax, N, timestep)
  for (m in 1:length(reductionStrategies)){
    percReduce <- reductionStrategies[m]
    parameters <- c(EFD, pEA, MDR, K_thr, a_reduceBiteRate, pMI, mu_th, PDR, b, timestep=timestep)
    state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
    out <- ode(y = state, times = times, func = seiseir_model_thr_reduceBiteRate, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$Site <- sites[l]
    out2$Intervention_strategy <- percReduce
    write.table(out2, file=traitFileName, row.names = F, sep = ",", col.names = !file.exists(traitFileName), append = T)
    cat("finished running ode for", sites[l], "intervention", m, "\n")
  }
}
