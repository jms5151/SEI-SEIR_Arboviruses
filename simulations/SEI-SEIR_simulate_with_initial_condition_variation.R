# Test initial conditions for SEI-SEIR model ----------------------------------------------------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(dplyr)
library(deSolve)

# load model 
source("models/SEI-SEIR_model_THR.R")

# load data 
source("models/SEI-SEIR_simulation_setup.R")

# run simulations using different initial conditions for state variables and save results -------------
initDF <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(initDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "Site", "InitCond")
fileName <- "output/Initital_condition_simulations.csv"
write.csv(initDF, fileName, row.names=F)

# run model simulations
for (l in 1:length(sites)){
  climateData2 <- subset(climateData, Site == sites[l])
  climateData2 <- climateData2[order(climateData2$Date),]
  climateData2 <- climateData2[complete.cases(climateData2),]
  temp <- climateData2$Temperature
  hum <- climateData2$SVPD
  rain <- climateData2$Two_week_rainfall
  Rmax <- 123
  if (unique(climateData2$country) == "Ecuador"){
    K_trh <- K_trh_quadratic
  } else {
    K_trh <- K_trh_briere
  }
  Date <- climateData2$Date
  H0 <- population[l]
  city <- sites[l]
  BR <- BRs[l]
  DR <- DRs[l]
  times <- seq(1,length(Date), by=1)
  M0 <- K_thr(temp[1], mean(rain), Rmax, H0, timestep)
  parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, timestep=timestep)
  for (k in 1:nrow(init.cond)){
    state <- c(M1 = init.cond$m1[k]*M0, M2 = init.cond$m2[k]*M0, M3 = init.cond$m3[k]*M0, S = init.cond$s[k]*H0, E = init.cond$e[k]*H0, I = init.cond$i[k]*H0, R = init.cond$r[k]*H0)
    out <- ode(y = state, times = times, func = seiseir_model_thr, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$Site <- sites[l]
    out2$IC <- c(init.cond$IC[k])
    write.table(out2, file=fileName, row.names = F, sep = ",", col.names = !file.exists(fileName), append = T)
    cat("finished running ode for", sites[l], "Initial conditions =", init.cond$IC[k], "\n")
  }
}
