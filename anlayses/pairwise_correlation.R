# calculate cross-correlation values and p-values between pairs of predictions and observations
rm(list=ls()) #remove previous variable assignments

# load library
library(tidyverse)

# load data
load("data/Cases.RData")
load("data/vector_data.RData")
mods <- read.csv("output/SEI-SEIR_simulations_THR.csv", head=T, stringsAsFactors = F)

# merge data
arboviruses <- merge(cases, mods, by=c("Site", "Date"), all.x = T)
aedes <- merge(vectors, mods, by=c("Site", "Date"), all.x = T)

# calculate total modeled mosquitoes
aedes$Mtot <- aedes$M1 + aedes$M2 + aedes$M3

# create empty data frame
df <- data.frame()

# calculate cross-correlation coefficients --------------------------
sites <- unique(mods$Site)
rain_fn <- unique(mods$Rain_function)
Ys <- c("cases", "vectors")

for (i in 1:length(sites)){
  for (j in 1:length(rain_fn)){
    for (k in 1:length(Ys)){
      if (Ys[k] == "cases"){
        pred <- "I"
        obs <- "arboviruses_positive"
        x <- arboviruses
      } else {
        pred <- "Mtot"
        obs <- "aedes_total"
        x <- aedes
      }
      # subset data
      x <- subset(x[,c("Site", "Rain_function", obs, pred)], Site == sites[i] & Rain_function == rain_fn[j])
      x <- x[complete.cases(x),]
      # calculate cross-correlation values for all possible lags
      ccfvalues <- ccf(x[,pred], x[,obs])
      # save data
      tmpdf <- data.frame("lag" = ccfvalues$lag, "CCF" = round(ccfvalues$acf, 2))
      tmpdf$Site = sites[i] 
      tmpdf$Y = Ys[k] 
      tmpdf$Rain_function = rain_fn[j]
      df <- rbind(df, tmpdf)
    }
  }
}

# subset to max CCF with lag for each site and response variable
df_max_ccf <- df %>%
  # remove lags beyond 4 months
  filter(lag >= -4 & lag <= 1) %>%
  group_by(Site, Y) %>%
  # keep site, y, and rain function combination with highest cross-correlation
  filter(CCF == max(CCF)) %>%
  group_by(Site, Y) %>%
  # if there are multiple lags with the same cross-correlation, keep row with shortest lag time
  filter(lag == max(lag))

# since all rain functions perform equally well for cases in Zaruma, keep rain function for cases in Zaruma that performs best for mosquitoes
zaruma_rain_fn_vectors <- df_max_ccf$Rain_function[df_max_ccf$Site=="Zaruma" & df_max_ccf$Y=="vectors"]
df_max_ccf <- df_max_ccf[-which(df_max_ccf$Site=="Zaruma" & df_max_ccf$Y=="cases" & df_max_ccf$Rain_function != zaruma_rain_fn_vectors),]

# modified chelton method -------------------------------------------
df_max_ccf$pvalue <- NA

for (l in 1:nrow(df_max_ccf)){
  if (df_max_ccf$Y[l] == "cases"){
    pred <- "I"
    obs <- "arboviruses_positive"
    x <- arboviruses
  } else {
    pred <- "Mtot"
    obs <- "aedes_total"
    x <- aedes
  }
  # subset data
  x <- subset(x[,c("Site", "Rain_function", obs, pred)], Site == df_max_ccf$Site[l] & Rain_function == df_max_ccf$Rain_function[l])
  x <- x[complete.cases(x),]
  # sample size of pairwise observations!
  N <- nrow(x) 
  # number of lags to use
  J <- round(N/5)
  # time series for comparison
  X <- x[,pred]
  Y <- x[,obs]
  # means of each time series
  Xbar <- mean(X)
  Ybar <- mean(Y)
  # calculate denominators for  eq. 6 (Pyper and Peterman 1998)
  ssX <- ssY <- NA
  for (ii in 1:length(Y)){
    ssX[ii] <- (X[ii]-Xbar)^2
    ssY[ii] <- (Y[ii]-Ybar)^2
  }
  denomX <- sum(ssX, na.rm=T)
  denomY <- sum(ssY, na.rm=T)
  # vectors of autocorrelation parameters
  rXX <- rYY <- NA 
  for (j in 1:J){
    ssX <- ssY <- NA
    Xl <- lag(X, j)
    Yl <- lag(Y, j)
    for(k in 1:(length(Y)-J)){ # loop through the numerator for eq. 6 (Pyper and Peterman 1998)
      ssX[k] <- (X[k]-Xbar)*(Xl[k]-Xbar)
      ssY[k] <- (Y[k]-Ybar)*(Yl[k]-Ybar)	
    }	
    numerX <- sum(ssX, na.rm=T) # numerator from eq. 6 (Pyper and Peterman 1998)
    numerY <- sum(ssY, na.rm=T)
    adj <- N/(N-j) # adjustor from eq. 7 (Pyper and Peterman 1998)
    rXX[j] <- adj*numerX/denomX
    rYY[j] <- adj*numerY/denomY
  }
  # calculate effective sample size (Nef)
  cross.prod <- NA
  for (j in 1:J){
    cross.prod[j] <- ((N-j)/N)*rXX[j]*rYY[j]
  }
  Nef <- 1/((1/N)+(2/N)*sum(cross.prod))
  if(Nef > N) Nef <- N # constrain to be no greater than N!
  r.ac <- df_max_ccf$CCF[l] # cor(X,Y, use="pair")
  t <- sqrt((Nef*r.ac^2)/(1-r.ac^2))
  p <- 2*pt(-abs(t), Nef-2)
  df_max_ccf$pvalue[l] <- round(p, 2)
}

# save data
save(df_max_ccf, file="output/CCF.RData")
