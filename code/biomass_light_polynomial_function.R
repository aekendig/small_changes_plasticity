#### set-up ####
library(deSolve)

#### biomass-light parameters ####

max_light <- 1800   # max light value (Drew's figure)
b <- 2   # max biomass size in grams (Drew's figure)


#### biomass-light functions ####

# linear function
bl_lin <- function(light, bio){
  
  # parameters
  a1 <- bio / max_light   # linear coefficient (Drew's figure)
  
  # equation
  biomass1 <- light * a1
  
  # return
  return(biomass1)
  
}

# quadratic function
bl_quad <- function(light, h, j, bio){
  
  # parameters
  x <- c(0, h, max_light)
  y <- c(0, j, bio)
  coeffs <- solve(cbind(1, x, x^2), y)
  c2 <- coeffs[1]
  b2 <- coeffs[2]
  a2 <- coeffs[3]
  
  # equation
  biomass2 <- a2 * light^2 + b2 * light + c2
  
  # return
  return(biomass2)
  
}

# main function
bl_fun <- function(light, poly, h = 0, j = 0, bio = 0){
  
  # choose function
  ifelse(poly == 1, return(bl_lin(light, bio)), return(bl_quad(light, h, j, bio)))
  
}
