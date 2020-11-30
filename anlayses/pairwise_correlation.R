# return maximum cross-correlation coefficient and lag 
cross_correlation_and_lags <- function(df, pred, obs){
  # calculate cross-correlation values for all possible lags
  ccfvalues <- ccf(df[, pred], df[, obs])
  # create data frame of output 
  tmpdf <- data.frame("lag" = ccfvalues$lag, "CCF" = round(ccfvalues$acf, 2))
  tmpdf <- tmpdf %>% filter(lag >= -4 & lag <= 0) %>% filter(CCF == max(CCF))
  return(list(tmpdf[,"CCF"], tmpdf[, "lag"]))
} 

# calculate p-value for cross-correlation using the modified chelton method
modified_chelton_method <- function(df, pred, obs, cross_correlation_val){
  # sample size of pairwise observations!
  N <- nrow(df) 
  # number of lags to use
  J <- round(N/5)
  # time series for comparison
  X <- df[, pred]
  Y <- df[, obs]
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
  if(Nef > N){
    Nef <- N # constrain to be no greater than N!
  } 
  t <- sqrt((Nef*r.ac^2)/(1-cross_correlation_val^2))
  p <- 2*pt(-abs(t), Nef-2)
  return(round(p, 2))
}