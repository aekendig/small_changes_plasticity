#### setup ####

# clear environment
rm(list=ls())

# load packages
library(tidyverse)

# figure settings
fig_theme <- theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.box.margin = margin(-10, -10, -10, -10))


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

# biomass-light parameters
max_light <- 1800   # max light value (Drew's figure)
b <- 2   # max biomass size in grams (Drew's figure)


#### biomass-light functions ####

# linear function
bl_lin <- function(light){
  
  a <- b / max_light   # linear coefficient (Drew's figure)
  
  biomass <- light * a
  
  return(biomass)
  
}

# positive function
bl_pos <- function(light){
  
  h <- max_light/2
  j <- b
  a <- -j/(h^2)
  biomass <- a * (light - h)^2 + j
  
  return(biomass)
  
}

# negative function
bl_neg <- function(light){
  
  h <- max_light/2
  j <- 0
  a <- (b-j)/(h^2)
  biomass <- a * (light - h)^2 + j
  
  return(biomass)
  
}

# dataframe of test functions
bl_test <- tibble(light_val = rep(seq(0, max_light, length.out = 200), 3),
                  fun_shape = rep(c("linear", "positive", "negative"), each = 200)) %>%
  mutate(biomass = case_when(fun_shape == "linear" ~ bl_lin(light_val),
                             fun_shape == "positive" ~ bl_pos(light_val),
                             fun_shape == "negative" ~ bl_neg(light_val)),
         fun_shape = fct_relevel(fun_shape, "linear", "positive"))

# figure of function test
pdf("output/model_v3/light_biomass_function_shapes.pdf")
ggplot(bl_test, aes(light_val, biomass, color = fun_shape)) +
  geom_line() +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 250, linetype = "dashed", color = "black") +
  xlab(expression(paste("Light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Biomass (g)") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme +
  theme(legend.position = c(0.8, 0.2))
dev.off()


#### intersection points ####

# make wide by function shape
bl_test_w <- pivot_wider(bl_test, names_from = fun_shape, values_from = biomass) %>%
  mutate(poslin = positive - linear,
         linneg = linear - negative)

# light values where functions intersect
light_pos <- filter(bl_test_w, poslin > 0) %>%
  filter(light_val == max(light_val))

light_neg <- filter(bl_test_w, linneg > 0) %>%
  filter(light_val == min(light_val))


#### simulation function ####

sim_fun <- function(light, fun_shape){
  
  # juvenile biomass as a function of light
  y_J <- ifelse(fun_shape == "linear", bl_lin(light), ifelse(fun_shape == "positive", bl_pos(light), bl_neg(light))) * 3 # assume this was over one month and increase for longer growing season
  
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
    mutate(biomass = juveniles * y_J + adults * y_A,
           trees = juveniles + adults)
  
  # return
  return(dfN)
}


#### test a few scenarios ####

# light is at 10 and biomass growth is a linear function of light
test_10_lin <- sim_fun(10, "linear") %>%
  as_tibble() %>%
  mutate(Light = 10,
         fun_shape = "linear")

# light is at 10 and biomass growth is a positive function of light
test_10_pos <- sim_fun(10, "positive") %>%
  as_tibble() %>%
  mutate(Light = 10,
         fun_shape = "positive")

# light is at 10 and biomass growth is an negative function of light
test_10_neg <- sim_fun(10, "negative") %>%
  as_tibble() %>%
  mutate(Light = 10,
         fun_shape = "negative")

# light is at 250 and biomass growth is a linear function of light
test_250_lin <- sim_fun(250, "linear") %>%
  as_tibble() %>%
  mutate(Light = 250,
         fun_shape = "linear")

# light is at 250 and biomass growth is a positive function of light
test_250_pos <- sim_fun(250, "positive") %>%
  as_tibble() %>%
  mutate(Light = 250,
         fun_shape = "positive")

# light is at 250 and biomass growth is a negative function of light
test_250_neg <- sim_fun(250, "negative") %>%
  as_tibble() %>%
  mutate(Light = 250,
         fun_shape = "negative")

# combine them
test_scen <- rbind(test_10_lin, test_10_pos, test_10_neg, test_250_lin, test_250_pos, test_250_neg) %>%
  mutate(fun_shape = fct_relevel(fun_shape, "linear", "positive"),
         Light = as.factor(Light))

# figures
pdf("output/model_v3/test_scenarios_time_series.pdf")
ggplot(test_scen, aes(time, seeds, color = fun_shape, linetype = Light)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Seeds") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme

ggplot(test_scen, aes(time, juveniles, color = fun_shape, linetype = Light)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Juveniles") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme

ggplot(test_scen, aes(time, adults, color = fun_shape, linetype = Light)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Adults") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme

ggplot(test_scen, aes(time, trees, color = fun_shape, linetype = Light)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Trees") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme

ggplot(test_scen, aes(time, biomass, color = fun_shape, linetype = Light)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Biomass (g)") +
  guides(color = guide_legend(title = "Function shape")) +
  fig_theme
dev.off()


#### final by resource change ####

# number of resource datapoints
res_points <- 200

# input dataframe
res_curve <- tibble(light = rep(seq(0, max_light, length.out = res_points), 3),
                    fun_shape = rep(c("linear", "positive", "negative"), each = res_points)) %>%
  mutate(seeds = NA,
         juveniles = NA,
         adults = NA,
         biomass = NA)

# final density for each input combination
for(i in 1:nrow(res_curve)){
  
  # run model
  temp_mod <- sim_fun(light = res_curve$light[i], fun_shape = res_curve$fun_shape[i])
  
  # add final timepoint to dataframe
  res_curve$seeds[i] <- temp_mod$seeds[simtime]
  res_curve$juveniles[i] <- temp_mod$juveniles[simtime]
  res_curve$adults[i] <- temp_mod$adults[simtime]
  res_curve$biomass[i] <- temp_mod$biomass[simtime]
  
}

# check 0 light values
filter(res_curve, light == 0)
# all 0

# modify dataframe
res_curve2 <- res_curve %>%
  mutate(trees = juveniles + adults)

# figure
pdf("output/model_v3/final_prediction_light_change.pdf")
ggplot(res_curve2, aes(light, trees, color = fun_shape)) +
  geom_line() +
  guides(color = guide_legend(title = "Function shape")) +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Final number of trees") +
  fig_theme

ggplot(res_curve2, aes(light, biomass, color = fun_shape)) +
  geom_line() +
  guides(color = guide_legend(title = "Function shape")) +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Final total biomass (g)") +
  fig_theme
dev.off()


#### difference from linear ####

# make wide by function shape
res_curve_w <- res_curve2 %>%
  select(-c(seeds:adults)) %>%
  pivot_wider(names_from = fun_shape,
              values_from = c(biomass, trees),
              names_glue = "{fun_shape}_{.value}") %>%
  mutate(linear_biomass_ch = linear_biomass - linear_biomass,
         linear_trees_ch = linear_trees - linear_trees,
         positive_biomass_ch = linear_biomass - positive_biomass,
         positive_trees_ch = linear_trees - positive_trees,
         negative_biomass_ch = linear_biomass - negative_biomass,
         negative_trees_ch = linear_trees - negative_trees)

# make long by change type
res_curve_l <- res_curve_w %>%
  select(light, linear_biomass_ch, positive_biomass_ch, negative_biomass_ch) %>%
  pivot_longer(cols = c(linear_biomass_ch, positive_biomass_ch, negative_biomass_ch),
               names_to = "fun_shape",
               values_to = "biomass_change") %>%
  mutate(fun_shape = sub("\\_.*", "", fun_shape)) %>%
  full_join(res_curve_w %>%
              select(light, linear_trees_ch, positive_trees_ch, negative_trees_ch) %>%
              pivot_longer(cols = c(linear_trees_ch, positive_trees_ch, negative_trees_ch),
                           names_to = "fun_shape",
                           values_to = "trees_change") %>%
              mutate(fun_shape = sub("\\_.*", "", fun_shape)))

# figure
pdf("output/model_v3/difference_linear_light_change.pdf")
ggplot(res_curve_l, aes(light, trees_change, color = fun_shape)) +
  geom_line() +
  guides(color = guide_legend(title = "Actual biomass-light relationship")) +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Tree difference due to linear approximation") +
  fig_theme +
  theme(legend.position = "bottom")

ggplot(res_curve_l, aes(light, biomass_change, color = fun_shape)) +
  geom_line() +
  guides(color = guide_legend(title = "Actual biomass-light relationship")) +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Biomass difference due to linear approximation (g)") +
  fig_theme +
  theme(legend.position = "bottom")
dev.off()


#### differences across resource levels ####

# divide data by function shape
res_curv_lin <- filter(res_curve2, fun_shape == "linear")
res_curv_pos <- filter(res_curve2, fun_shape == "positive")
res_curv_neg <- filter(res_curve2, fun_shape == "negative")

# compute the differences
res_diff_lin <- tibble(light_diff = as.numeric(dist(res_curv_lin$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(res_curv_lin$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(res_curv_lin$trees, method = "manhattan")),
                       fun_shape = "linear")
res_diff_pos <- tibble(light_diff = as.numeric(dist(res_curv_pos$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(res_curv_pos$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(res_curv_pos$trees, method = "manhattan")),
                       fun_shape = "positive")
res_diff_neg <- tibble(light_diff = as.numeric(dist(res_curv_neg$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(res_curv_neg$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(res_curv_neg$trees, method = "manhattan")),
                       fun_shape = "negative")

# combine
# divide trait differences be environmental differences
res_diff <- full_join(res_diff_lin, res_diff_pos) %>%
  full_join(res_diff_neg) %>%
  mutate(biomass_espi = biomass_diff / light_diff,
         trees_espi = trees_diff / light_diff,
         light_interval = cut_interval(light_diff, n = 100))

# light diff values
length(unique(res_diff$light_diff))

# figure
pdf("output/model_v3/espi.pdf")
ggplot(res_diff, aes(light_diff, biomass_espi)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ fun_shape, scales = "free") +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass ESPI (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme

ggplot(res_diff, aes(light_interval, biomass_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass ESPI (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme +
  theme(axis.text.x = element_blank())

ggplot(res_diff, aes(light_diff, trees_espi)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ fun_shape, scales = "free") +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Trees ESPI (number (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme

ggplot(res_diff, aes(light_interval, trees_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  fig_theme +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Trees ESPI (number (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  theme(axis.text.x = element_blank())
dev.off()


#### ESPI overlapping ranges ####

# divide data by function shape
ov_res_curv_lin_pos <- filter(res_curve2, fun_shape == "linear" & light <= light_pos$light_val)
ov_res_curv_lin_neg <- filter(res_curve2, fun_shape == "linear" & light >= light_neg$light_val)
ov_res_curv_pos <- filter(res_curve2, fun_shape == "positive" & light <= light_pos$light_val)
ov_res_curv_neg <- filter(res_curve2, fun_shape == "negative" & light >= light_neg$light_val)

# compute the differences
ov_res_diff_lin_pos <- tibble(light_diff = as.numeric(dist(ov_res_curv_lin_pos$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(ov_res_curv_lin_pos$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(ov_res_curv_lin_pos$trees, method = "manhattan")),
                       fun_shape = "linear")
ov_res_diff_lin_neg <- tibble(light_diff = as.numeric(dist(ov_res_curv_lin_neg$light, method = "manhattan")),
                              biomass_diff = as.numeric(dist(ov_res_curv_lin_neg$biomass, method = "manhattan")),
                              trees_diff = as.numeric(dist(ov_res_curv_lin_neg$trees, method = "manhattan")),
                              fun_shape = "linear")
ov_res_diff_pos <- tibble(light_diff = as.numeric(dist(ov_res_curv_pos$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(ov_res_curv_pos$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(ov_res_curv_pos$trees, method = "manhattan")),
                       fun_shape = "positive")
ov_res_diff_neg <- tibble(light_diff = as.numeric(dist(ov_res_curv_neg$light, method = "manhattan")),
                       biomass_diff = as.numeric(dist(ov_res_curv_neg$biomass, method = "manhattan")),
                       trees_diff = as.numeric(dist(ov_res_curv_neg$trees, method = "manhattan")),
                       fun_shape = "negative")

# combine
# divide trait differences be environmental differences
ov_res_diff_pos <- full_join(ov_res_diff_lin_pos, ov_res_diff_pos) %>%
  mutate(biomass_espi = biomass_diff / light_diff,
         trees_espi = trees_diff / light_diff,
         light_interval = cut_interval(light_diff, n = 100))

ov_res_diff_neg <- full_join(ov_res_diff_lin_neg, ov_res_diff_neg) %>%
  mutate(biomass_espi = biomass_diff / light_diff,
         trees_espi = trees_diff / light_diff,
         light_interval = cut_interval(light_diff, n = 100))

# figure
pdf("output/model_v3/espi_overlap.pdf")
ggplot(ov_res_diff_pos, aes(light_interval, biomass_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass ESPI (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme +
  theme(axis.text.x = element_blank())

ggplot(ov_res_diff_pos, aes(light_interval, trees_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  fig_theme +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Trees ESPI (number (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  theme(axis.text.x = element_blank())

ggplot(ov_res_diff_neg, aes(light_interval, biomass_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass ESPI (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme +
  theme(axis.text.x = element_blank())

ggplot(ov_res_diff_neg, aes(light_interval, trees_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ fun_shape) +
  fig_theme +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Trees ESPI (number (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  theme(axis.text.x = element_blank())
dev.off()
