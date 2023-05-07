rm(list=ls())
library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(ggplot2)
library(rootSolve)

####################################################################################################################
# Calibration module
####################################################################################################################

# I. Totoaba growth #######
K = 20226 # Carrying capacity in tonnes
avg_weight = 28 # Individual average weight in jilograms
K_ind = K*1000/avg_weight # Carrying capacity in individuals
r = .23 # Population growth rate

# II. Demand characteristics ########
# Assume a linear inverse demand function p(q) = alpha - beta * q
# Where q : totoaba buche (unitless) of mass m in grams
#
alpha = 30000
beta = 50

# III. Cost & engineering characteristics ######
# A. For poachers : 
# Profit is : wage/fish * number of fish - W*Effort^2
W = 3588.436
# Catchability 
sigma = .1

# B. For traders : 
# c is cost of bribing, transport etc
c = 1208.372 


# Save parameters from calibration to data folder ####
parameters = data.frame(K, avg_weight, K_ind, r, alpha, beta, W, sigma, c)
write.csv(parameters, file = paste0(here(),'/data/params_simul.csv'))
