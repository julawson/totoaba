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
K = 20226 # Carrying capacity (in tonnes)
avg_weight = 28 # Individual average weight (in kilograms)
K_ind = K*1000/avg_weight # Carrying capacity in individuals (unitless)
r = 0.23 # Population growth rate (unitless) -- the value I found in the Upsides database is 0.039 (much lower)

# II. Demand characteristics ########
# Assume a linear inverse demand function p(q) = alpha - beta * q
# Where q : totoaba buche (unitless) of mass m in grams
#
alpha = 30000 #choke price (in USD per individual)
beta = 50 

# III. Cost & engineering characteristics ######
# A. For poachers : 
# Profit is : wage/fish * number of fish - W*Effort^2
W = 3588.436 # cost parameter (USD per unit of fish)
# Catchability 
sigma = .1 #(unitless)

# B. For traders : 
# c is cost of bribing, transport etc
c = 1208.372 # (in USD)


# Save parameters from calibration to data folder ####
parameters = as.data.frame(K, avg_weight, K_ind, r, alpha, beta, W, sigma, c)
write_csv(parameters, here("data","params_simul.csv"))
