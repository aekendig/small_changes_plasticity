#### setup ####

# clear environment
rm(list=ls())

# load packages
# library(plyr)
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

# source functions
source("code/seeds_juvenile_adult_discrete_model_simulation.R")
source("code/biomass_light_polynomial_function.R")


#### biomass-light relationships ####

# light values
light_val = seq(10, max_light, by = 10)
light_n <- length(light_val)
                
# linear relationships
bl_lin_hi <- bl_fun(light = light_val, poly = 1, bio = b)
bl_lin_lo <- bl_fun(light = light_val, poly = 1, bio = b/2)

# quadratic relationships
bl_up_hi <- bl_fun(light = light_val, poly = 2, h = max_light/2, j = 3, bio = b)
bl_down_hi <- bl_fun(light = light_val, poly = 2, h = max_light/2, j = 2, bio = b)
bl_up_lo <- bl_fun(light = light_val, poly = 2, h = max_light/2, j = 2, bio = b/2)
bl_down_lo <- bl_fun(light = light_val, poly = 2, h = max_light/2, j = 1, bio = b/2)

# combine
bl_curves <- tibble(light = rep(light_val, 6),
                    biomass = c(bl_lin_hi, bl_lin_lo, bl_up_hi, bl_down_hi, bl_up_lo, bl_down_lo),
                    approach = rep(c(rep("traditional", 2), rep("multiple points", 4)), each = light_n),
                    stress_change = rep(c("none", "lower slope", "none", "lower peak", "lower slope", "lower slope\nand peak"), each = light_n)) %>%
  mutate(approach = fct_relevel(approach, "traditional"),
         stress_change = fct_relevel(stress_change, "none", "lower slope", "lower peak"))

# points
bl_points <- bl_curves %>%
  filter((approach == "traditional" & light %in% c(10, max_light)) | (approach == "multiple points" & light == 500)) %>%
  mutate(app_stress = paste(approach, stress_change, sep = " ") %>%
           fct_relevel("traditional none", "traditional lower slope"))

# figure
pdf("output/model_v4/light_biomass_function_shapes.pdf")
ggplot(bl_curves, aes(light, biomass, color = stress_change)) +
  geom_line(aes(linetype = approach)) +
  geom_point(data = bl_points, size = 3, aes(shape = stress_change, color = stress_change, fill = app_stress)) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Approach") +
  scale_shape_manual(values = c(21, 21, 22, 22), name = "Stress effect") +
  scale_color_manual(values = c("#286C7B", "#A74A44", "#A2CEDE", "#E93F33"), name = "Stress effect") +
  scale_fill_manual(values = c("#286C7B", "#A74A44", rep("white", 4)), guide = F) +
  xlab(expression(paste("Light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Biomass (g)") +
  fig_theme +
  theme(legend.position = c(0.16, 0.8))
dev.off()


#### ESPI ####

# ESPI function
espi_fun <- function(dat){
  
  # compute differences
  dat_out <- tibble(light_diff = as.numeric(dist(dat$light, method = "manhattan")),
                    biomass_diff = as.numeric(dist(dat$biomass, method = "manhattan")),
                    approach = unique(dat$approach),
                    stress_change = unique(dat$stress_change)) %>%
    mutate(biomass_espi = biomass_diff / light_diff)
  
  # output
  return(dat_out)
  
}

# apply function
bl_espi <- plyr::ddply(bl_curves, c("approach", "stress_change"), espi_fun) %>%
  as_tibble() %>%
  mutate(light_interval = cut_interval(light_diff, n = 100))

# figure
pdf("output/model_v4/biomass_espi_separate.pdf")
ggplot(bl_espi, aes(light_diff, biomass_espi)) +
  geom_point(alpha = 0.6) +
  facet_grid(stress_change ~ approach) +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass plasticity (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme

ggplot(bl_espi, aes(light_interval, biomass_espi)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_grid(stress_change ~ approach) +
  xlab(expression(paste("Binned change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass plasticity (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme +
  theme(axis.text.x = element_blank())
dev.off()

pdf("output/model_v4/biomass_espi.pdf")
ggplot(bl_espi, aes(light_interval, biomass_espi, color = stress_change, shape = approach)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  scale_color_manual(values = c("#286C7B", "#A74A44", "#A2CEDE", "#E93F33"), name = "Stress effect") +
  scale_shape_manual(values = c(19, 1), name = "Approach") +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab(expression(paste("Biomass plasticity (g (", mu, "moles ", m^-2, " ", sec^-1, ")"^-1, ")", sep = ""))) +
  fig_theme +
  theme(axis.text.x = element_blank(),
        legend.position = c(0.8, 0.8))
dev.off()


#### final by resource change ####

# input dataframe
res_curve <- bl_curves %>%
  mutate(seeds = NA,
         juveniles = NA,
         adults = NA,
         total_biomass = NA)

# final density for each input combination
for(i in 1:nrow(res_curve)){
  
  # run model
  temp_mod <- sim_fun(bio = res_curve$biomass[i])
  
  # add final timepoint to dataframe
  res_curve$seeds[i] <- temp_mod$seeds[simtime]
  res_curve$juveniles[i] <- temp_mod$juveniles[simtime]
  res_curve$adults[i] <- temp_mod$adults[simtime]
  res_curve$total_biomass[i] <- temp_mod$total_biomass[simtime]
  
}

# modify dataframe
res_curve2 <- res_curve %>%
  mutate(trees = juveniles + adults)

# density function
trees_fun <- function(dat){
  
  # compute differences
  dat_out <- tibble(light_diff = as.numeric(dist(dat$light, method = "manhattan")),
                    trees_diff = as.numeric(dist(dat$trees)),
                    approach = unique(dat$approach),
                    stress_change = unique(dat$stress_change))
  
  # output
  return(dat_out)
  
}

# apply function
tree_change <- plyr::ddply(res_curve2, c("approach", "stress_change"), trees_fun) %>%
  as_tibble() %>%
  mutate(light_interval = cut_interval(light_diff, n = 100))

pdf("output/model_v4/trees_change.pdf")
ggplot(res_curve2, aes(light, log(trees), color = stress_change, linetype = approach)) +
  geom_line() +
  scale_color_manual(values = c("#286C7B", "#A74A44", "#A2CEDE", "#E93F33"), name = "Stress effect") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Approach") +
  xlab(expression(paste("Light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Log-transformed final population size") +
  fig_theme +
  theme(legend.position = c(0.8, 0.4))

ggplot(tree_change, aes(light_interval, trees_diff, color = stress_change, shape = approach)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  scale_color_manual(values = c("#286C7B", "#A74A44", "#A2CEDE", "#E93F33"), name = "Stress effect") +
  scale_shape_manual(values = c(19, 1), name = "Approach") +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Change in population size") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        legend.position = c(0.16, 0.8))

ggplot(tree_change, aes(light_interval, trees_diff, color = stress_change, shape = approach)) +
  stat_summary(geom = "point", fun = "mean") +
  scale_color_manual(values = c("#286C7B", "#A74A44", "#A2CEDE", "#E93F33"), name = "Stress effect") +
  scale_shape_manual(values = c(19, 1), name = "Approach") +
  xlab(expression(paste("Change in light availability (", mu, "moles ", m^-2, " ", sec^-1, ")", sep = ""))) +
  ylab("Change in population size") +
  fig_theme +
  theme(axis.text.x = element_blank(),
        legend.position = c(0.16, 0.8))
dev.off()

