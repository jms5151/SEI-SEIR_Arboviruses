# SEI-SEIR model with temperature, humidity, and rainfall --------------------------------
seiseir_model_thr <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    dM1 <- EFD(temp[t])*pEA(temp[t])*MDR(temp[t])*mu_th(temp[t], hum[t], timestep)^(-1)*(M1+M2+M3)*max((1-((M1+M2+M3)/K_thr(temp[t], rain[t], Rmax, (S+E+I+R), timestep))),0)-(a(temp[t])*pMI(temp[t])*(I/(S+E+I+R))+mu_th(temp[t], hum[t], timestep))*M1
    dM2 <- a(temp[t])*pMI(temp[t])*(I/(S+E+I+R))*M1-(PDR(temp[t])+mu_th(temp[t], hum[t], timestep))*M2
    dM3 <- PDR(temp[t])*M2-mu_th(temp[t], hum[t], timestep)*M3
    dS <- -a(temp[t])*b(temp[t])*(M3/(S+E+I+R))*S + BR*(S/1000)/360 - DR*(S/1000)/360 + ie*(S+E+I+R) - ie*S
    dE <- a(temp[t])*b(temp[t])*(M3/(S+E+I+R))*S-(1.0/5.9)*E - DR*(E/1000)/360 - ie*E
    dI <- (1.0/5.9)*E-(1.0/5.0)*I - DR*(I/1000)/360 - ie*I
    dR <- (1.0/5.0)*I - DR*(R/1000)/360 - ie*R
    list(c(dM1, dM2, dM3, dS, dE, dI, dR))
  })
}    

# This is the general function for the Briere fit.
briere <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*x*(x-T0)*sqrt(Tm-x)
}

# This is the general function for the quadratic fit. 
quadratic <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*(x-T0)*(x-Tm)
}

# This is the general function for the inverted quadratic fit.
inverted_quadratic <- function(x, c, T0, Tm, timestep){
  if((x < T0) | (x > Tm))
    1.0/timestep
  else
    1.0/(c*(x-T0)*(x-Tm))
}

# eggs per female per day
EFD <- function(temp){
  briere(temp, EFD_c, EFD_T0, EFD_Tm)
}

# probability egg to adult survival
pEA <- function(temp){
  quadratic(temp, pEA_c, pEA_T0, pEA_Tm)
}

# mosquito development rate (1/larval development period)
MDR <- function(temp){
  briere(temp, MDR_c, MDR_T0, MDR_Tm)
}

# biting rate
a <- function(temp){
  briere(temp, a_c, a_T0, a_Tm)
}

# probability	of mosquito	infection per	bite	on	an	infectious	host
pMI <- function(temp){
  briere(temp, pMI_c, pMI_T0, pMI_Tm)
}

# adult mosquito mortality rate (1/adult lifespan)
mu_th <- function(temp, hum, timestep){
  if (hum <= 1){
    inverted_quadratic(temp, mu_th_c, mu_th_T0, mu_th_Tm, timestep)+(1-(0.01256 + 2.00893*hum))*0.005
  } else {
    inverted_quadratic(temp, mu_th_c, mu_th_T0, mu_th_Tm, timestep)+(1-(1.2248 + 0.2673*hum))*0.01
  }
}

# parasite development rate
PDR <- function(temp){
  briere(temp, PDR_c, PDR_T0, PDR_Tm)
}

# transmission competence: probability of human	infection	per	bite	by	an	infectious mosquito
b <- function(temp){
  briere(temp, b_c, b_T0, b_Tm)
}

# carrying capacity
carrying_capacity_th <- function(temp, h0, T0, EA, N, timestep){
  kappa <- 8.617e-05; # Boltzmann constant 
  alpha <- (EFD(T0)*pEA(T0)*MDR(T0)*mu_th(T0, h0, timestep)^(-1)-mu_th(T0, h0, timestep))/(EFD(T0)*pEA(T0)*MDR(T0)*mu_th(T0, h0, timestep)^(-1))
  (alpha*N*exp(-EA*((temp-T0)^2)/(kappa*(temp+273.0)*(T0+273.0))))
}

# carrying capacity with temperature, humidity, and rainfall
K_thr_briere <- function(temp, rain, Rmax, N, timestep){
  R0 <- 1
  if((rain < R0) | (rain > Rmax)){
    max(0.01*carrying_capacity_th(temp, 100, 29.0, 0.05, N, timestep), 1000)
  }
  else {
    c <- 7.86e-05
    max(carrying_capacity_th(temp, 100, 29.0, 0.05, N, timestep)*c*rain*(rain-R0)*sqrt(Rmax-rain)*0.268 + 0.001, 1000)
  }
}

K_thr_quadratic <- function(temp, rain, Rmax, N, timestep){
  R0 <- 1
  if((rain < R0) | (rain > Rmax)){
    max(0.01*carrying_capacity_th(temp, 100, 29.0, 0.05, N, timestep), 1000)
  }
  else {
    c <- -5.99e-03
    max(carrying_capacity_th(temp, 100, 29.0, 0.05, N, timestep)*(c*(rain-R0)*(rain-Rmax))*0.045 + 0.001, 1000)
  }
}

K_thr_inverse <- function(temp, rain, Rmax, N, timestep){
  if (rain < 1){
    rain <- 1
  } 
  max(carrying_capacity_th(temp, 100, 29.0, 0.05, N, timestep)*(1/rain), 1000)
}
