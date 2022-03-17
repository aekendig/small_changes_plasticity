#### set-up ####

# load packages
library(tidyverse)
library(deSolve)


#### model ####

# define model
logistic_mod <- function(time, y, parms){
  
  with(as.list(c(y, parms)), {
    
    dN <- r * N * (1 - N/K)
    
    list(c(dN))
    
  })
  
}

# set values
y0 <- c(N = 0.1)
times <- seq(0, 100, 1)
parms_hi <- c(r = 0.1, K = 20)
parms_lo <- c(r = 0.1, K = 10)

# simulations
out_hi <- ode(y0, times, logistic_mod, parms_hi)
out_lo <- ode(y0, times, logistic_mod, parms_lo)

# quick plots
plot(out_hi)
lines(out_lo)


#### figure ####

# combine data
out <- out_hi %>%
  as_tibble() %>%
  mutate(Resource = "high",
         time = as.numeric(time),
         N = as.numeric(N)) %>%
  full_join(out_lo %>%
              as_tibble() %>%
              mutate(Resource = "low",
                     time = as.numeric(time),
                     N = as.numeric(N)))

# figure
fig <- ggplot(out, aes(x = time, y = N, color = Resource)) +
  geom_line()

ggsave("output/logistic_model_resources.pdf", fig, width = 5, height = 3)
