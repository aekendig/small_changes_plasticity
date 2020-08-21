#### equations ####

# seeds
# S[t+1] = s(1-g)S[t] + c(y_J*J[t] + y_A*A[t])/(1 + alpha(y_J*J[t] + y_A*A[t]))
# seeds produced this year = seeds that didn't germinate and survived in the seed bank + seed produced this year, which depends on the size of adults and juveniles and intraspecific competition

# juveniles
# J[t+1] = gS[t] + (1-m)J[t]
# juveniles this year = seeds that germinated and survived + juveniles that didn't mature

# adults
# A[t+1] = mJ[t] + kA[t]
# adults this year = juveniles that matured and adults that survived from last year


#### parameters ####

# initial conditions
S0 <- 100   # initial seed bank
J0 <- 0   # initial juvenile trees
A0 <- 0   # initial adult trees

# simulation conditions
simtime <- 500  # years of simulation

# plant traits (all per year)
s <- 0.8   # survival in the seed bank (fake value)
g <- 0.4   # proportion of seeds that germinate (fake value)
c <- 39000 / 1000   # seeds per gram (https://www.cabi.org/isc/datasheet/18119)
alpha <- 0.2   # reduction in seed production based on plant biomass (fake value)
m <- 0.8   # juvenile survival and maturation (fake value)
k <- 0.9   # adult survival (fake value)
q <- 4   # conversion from juvenile to adult biomass (fake value)


#### simulation function ####

sim_fun <- function(bio){
  
  # juvenile biomass as a function of light
  y_J <- bio * 3 # assume this was over one month and increase for longer growing season
  
  # adult biomass
  y_A <- q * y_J
  
  # initialize populations
  S <- rep(NA,simtime)
  J <- rep(NA,simtime)
  A <- rep(NA,simtime)
  
  S[1] <- S0
  J[1] <- J0
  A[1] <- A0
  
  # simulate population dynamics
  for(t in 1:(simtime - 1)){	
    
    # population size
    S[t+1] = s * (1 - g) * S[t] + c * (y_J * J[t] + y_A * A[t]) / (1 + alpha * (y_J * J[t] + y_A * A[t]))
    J[t+1] = g * S[t] + (1 - m) * J[t]
    A[t+1] = m * J[t] + k * A[t]
    
    # correct to prevent negative numbers
    S[t+1] = ifelse(S[t+1] < 1, 0, S[t+1])
    J[t+1] = ifelse(J[t+1] < 1, 0, J[t+1])
    A[t+1] = ifelse(A[t+1] < 1, 0, A[t+1])
  }
  
  # save data
  dfN = tibble(time = 1:simtime, seeds = S, juveniles = J, adults = A) %>%
    mutate(total_biomass = juveniles * y_J + adults * y_A,
           trees = juveniles + adults)
  
  # return
  return(dfN)
}