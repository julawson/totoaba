rm(list=ls())
library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(ggplot2)
library(rootSolve)

####################################################################################################################
# Simulation module
####################################################################################################################

# Step 1 : Initial simulation from calibration #####
parameters = read.csv(paste0(here(), '/data/params_simul.csv'))

x = seq(1, 
        parameters$K_ind,
        10)


# I. Growth ####
growth = function(x){
  r = parameters$r
  K = parameters$K_ind
  y = r*x*(1-x/K)
  return(y)
}

# II. SMonopoly ####
# A. Poachers wage
monopoly_wage = function(x){
  alpha = parameters$alpha
  c = parameters$c
  W = parameters$W
  beta = parameters$beta
  sigma = parameters$sigma
  y = 2*W*(alpha - c)/(2*sigma^2*beta*x^2 + 2*W)
  return(y)
}

monopoly_poaching = function(x){
  alpha = parameters$alpha
  c = parameters$c
  W = parameters$W
  
  y = (alpha - c)/(2*parameters$beta * parameters$sigma^2 * x^2 + 2*W)
  return(y)
}

results = data.frame(x)%>%
  mutate(growth = growth(x),
         monopoly = monopoly_poaching(x), 
         diff = growth - monopoly)

# Identify steady states : should replicate initial x
results %>% subset(diff < .1 & diff >-1.) 
# Would be nice to have something that extracts x or solves diff for steady state
# Output is stock, salary, quantities and profit
# III. Cournot competition ####
### A. Poachers wage ####
wage_cournot = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (2*W*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/((sigma^2)*(x^2)*(4*beta_f*beta_w - gamma^2) + 4*beta_f*W)
  return(y)
}

### B. Quantities ####
q_cournot_wild= function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/(4*beta_f*W + (sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)))
  return(y)
}

q_cournot_farmed= function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_w*(alpha_f - v) - gamma*(alpha_w-c)) + 2*W*(alpha_f-v))/((4*beta_f*beta_w - gamma^2)*(sigma^2)*(x^2)+ 4*W*beta_f)
  return(y)
}

### C. Solve for steady state with new diff : growth - q_cournot_wild ####

### D. Outputs ####
# Compute quantity, profits

# IV. Betrand Competition #####
# A. Poachers wage
wage_betrand = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = (2*W*b_w*( b_f*(2*a_w + e*v) + c*((e^2) - 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + 2*W*b_w*(2*b_w*b_f - (e^2)))
  return(y)
}

# B. Quantities
q_betrand_wild = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v) + c*((e^2)- 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + 2*W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

q_bertrand_farmed = function(s){
  y = b_f *(2*b_w*a_f + v*(e^2 - 2*b_w*b_f) + e*(a_w + (s+c)*b_w))/(4*b_f*b_w - e^2)
  return(y)
}
# C. Steady states 
# Outputs being stock

# D. Outputs




# Step 2 : Sensitivity analysis/scenario analysis
# Can factor in uncertainty about parameters
# think of expand.grid to test multiple combinations

# Loop over each set of parameters
for (i in length(params_df)){
  params = params_df[i,]
  # Run simulation 
  
  # Find a clever way to record steady states, i.e when diff gets to 0
}