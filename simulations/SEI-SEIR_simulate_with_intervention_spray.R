# SEI-SEIR model simulations with intervention: reduce mosquito population --------------
rm(list=ls()) #remove previous variable assignments

# load packages
library(deSolve)

# load data 
source("Codes/SEI-SEIR_model_with_intervention_spray.R")
source("Codes/SEI-SEIR_simulation_setup.R")

# create intervention strategy
timespan <- seq.Date(from=as.Date("2016-10-01"), to=as.Date("2018-03-31"), by="days")
InterveneDates <- seq.Date(from=as.Date("2017-01-01"), to=as.Date("2017-12-31"), by="months")

simInt <- data.frame(matrix(ncol=2, nrow=0)); colnames(simInt) <- c("Date", "Intervention")
reductionStrategies <- c(0.10, 0.50, 0.90)

for (i in 1:length(InterveneDates)){
  for (j in 1:length(reductionStrategies)){
    df <- data.frame(Date=timespan, Intervention=0)
    df$Intervention[df$Date==InterveneDates[i]] <- reductionStrategies[j]
    simInt <- rbind(simInt, df)
  }
}

simInt$simID <- sort(rep(1:(12*length(reductionStrategies)),length(timespan)))
simInt$simID <- simInt$simID+1
nullYear <- data.frame(Date=timespan, Intervention=0, simID=1) # add year with no intervention
simInt <- rbind(nullYear, simInt)

# run simulations
traitDF <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(traitDF) <- c("time", "M1", "M2", "M3", "S", "E", "I", "R", "Date", "Site", "SimID", "Intervention_strategy", "Intervention_date")
traitFileName <- "Concatenated_Data/model_simulations/SEI-SEIR_simulations_with_intervention_spray.csv"
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
  M0 <- K_thr(temp[1], hum[1], mean(rain), Rmax, N)
  for (m in 1:length(unique(simInt$simID))){
    Intervention <- simInt$Intervention[simInt$simID==m]
    parameters <- c(EFD, pEA, MDR, K_thr, a, pMI, mu_th, PDR, b, intv, timestep=timestep)
    state <- c(M1 = startIC$m1*M0, M2 = startIC$m2*M0, M3 = startIC$m3*M0, S = startIC$s*N, E = startIC$e*N, I = startIC$i*N, R = startIC$r*N)
    out <- ode(y = state, times = times, func = seiseir_model_thr_spray, parms = parameters, method="rk4", atol = 1e-14, rtol = 1e-14, hini = timestep)
    out2 <- as.data.frame(out)
    out2$Date <- Date
    out2$Site <- sites[l]
    out2$SimID <- m
    out2$Intervention_strategy <- ifelse(length(unique(Intervention))>0, max(Intervention), 0)
    out2$Intervention_date <- ifelse(unique(out2$Intervention_strategy)==0, "None", as.character(timespan[which(Intervention==max(Intervention))]))
    write.table(out2, file=traitFileName, row.names = F, sep = ",", col.names = !file.exists(traitFileName), append = T)
    cat("finished running ode for", sites[l], "intervention", m, "\n")
  }
}
