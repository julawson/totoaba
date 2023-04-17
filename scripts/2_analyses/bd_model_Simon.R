library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(rootSolve)

### I. Parameters #####
# Define the parameters for the biological resource
r <- 0.16 #intrinsic rate of population increase, estimated
k <- 100000 #carrying capacity, number of individuals
sigma_min <- 0.000167 #catchability, ranging from 0.000167 to 0.000256
sigma_max <- 0.000256
sigma_own = 0.0003677

# Define the parameters for the poaching effort
n <- 0.1
P <- 221 #exogenous market price per unit of output
W <- 5 #cost parameter, range is 2.18 to 7.6, midpoint 5
theta <- 2 #cost parameter, assumed to be 2
z <- 2 * W
s <- 25 #price paid to poachers, ranges from $15-25 
#Parameters for the trader's payoff function
alpha <- 6182 #parameter associated with the inverse demand function
beta <- 2.13 #parameter associated with the inverse demand function
c <- 50 #costs associated with the transportation and final sale of wild animal products 


# Initial conditions
initial_x <- 2600 #number of individuals
initial_E <- 0.000104 #presumably catch per unit effort, which is 0.000104 to 0.15 for rhinos

# Time points for the simulation
time_points <- seq(0, 50, 0.1)


### II. Monopoly set-up #####
x = seq(0,k,1) # Generate space for x

# Define growth function
growth <- function(x,r,k){
  return(r*x*(1-x/k))
}

#Initiate results dataframe
results = data.frame(x)

results$growth = growth(x,r,k)  

# Define harvest function in Monopoly
harvest_monop <- function(alpha, c, sigma, beta, z, x){
  return((alpha - c)*(sigma^2)*(x^2)/(2*beta*(sigma^2)*(x^2)+z))
}
# Set results 
results$harvest_monop_min = harvest_monop(alpha, c, sigma_min, beta, z, x)
results$harvest_monop_max = harvest_monop(alpha, c, sigma_max, beta, z, x)
results$harvest_monop_own = harvest_monop(alpha, c, sigma_own,beta,z,x)

# Compute growth - harvest
results = results %>% mutate(steady_monop_min = growth - harvest_monop_min,
                             steady_monop_max = growth - harvest_monop_max,
                             steady_monop_own = growth - harvest_monop_own)

# Find where growth = harvest to find steady states (sign change is necessary)
results %>% subset(steady_monop_min>-0.1 & steady_monop_min<0.1)
# Steady state is at 90,120 in this case

results %>% subset(steady_monop_max>-0.1 & steady_monop_max<0.1)
# Steady state is at 90,054

results %>% subset(steady_monop_own>-0.1 & steady_monop_own<0.1)
# 3 steady states : 2,625; 7,346; 90,029

# Values in the paper are : 
# 2,625
# 7,346
# 90,030


# Plot for visual exploration:
results %>% ggplot(aes(x)) + 
  geom_point(aes(y=growth, colour = 'growth')) + 
  geom_point(aes(y=harvest_monop_min, colour='harvest_monop_min')) +
  geom_point(aes(y=harvest_monop_max, colour='harvest_monop_max'))+
  geom_point(aes(y=harvest_monop_own, color = 'harvest_monop_sens'))

### III. Cournot model ####
# New parameters: 
v = 1000 # farming costs
gamma = 0.75 # substitutability
alpha_f = alpha # Demand intercept
alpha_w = alpha # 
beta_w = beta # Demand slope
beta_f = beta

# Reaction functions : 
r_f = function(q_w,alpha_f,gamma,w,v, beta_f){
  return((alpha_f - gamma*q_w - w)/2*beta_f)
}
r_w = function(q_f,alpha_w,gamma,s,c,beta_w){
  return((alpha_w - gamma*q_f - s - c)/2*beta_w)
}


# Optimal remuneration
delta = 4*beta_f*beta_w-gamma^2

s_c_min = (W*(2*beta_f*alpha_w-gamma*(alpha_f-v)-2*c*beta_f))/(sigma_min^2*x^2*delta+2*beta_f*W)
s_c_max = (W*(2*beta_f*alpha_w-gamma*(alpha_f-v)-2*c*beta_f))/(sigma_max^2*x^2*delta+2*beta_f*W)
s_c_own = (W*(2*beta_f*alpha_w-gamma*(alpha_f-v)-2*c*beta_f))/(sigma_own^2*x^2*delta+2*beta_f*W)


# Quantities
# Farming production in Cournot
q_c_F = function(x, sigma, alpha_w, beta_w, gamma, c, v, W, alpha_f, delta){
  y = (sigma^2*x^2*(2*alpha_f*beta_w + gamma*(alpha_w-c) - 2*beta_w*v) + W*(alpha_f-v))/(2*beta_f*W + sigma^2*x^2*delta)
  return(y)
}

# Harvesting of wild population in Cournot
q_c_w = function(x,sigma,alpha_w, beta_f, alpha_f, gamma, c, v, W, delta){
  y = (sigma^2)*(x^2)*(2*alpha_w*beta_f - gamma*(alpha_f-v) - 2*beta_f*c)/(2*beta_f*W + (sigma^2)*(x^2)*delta)
  return(y)
}

q_mod_c_w = function(beta_f,beta_w, alpha_w, s, c, gamma, v){
  y = (2*beta_f*(alpha_w-s-c)-gamma*(alpha_w-v))/(4*beta_w*beta_f - gamma^2)
  return(y)
}

results = results %>% mutate(harvest_cournot_min = q_c_w(x, sigma_min, alpha_w, beta_f, alpha_f, gamma, c, v, W, delta),
                             harvest_cournot_max = q_c_w(x, sigma_max, alpha_w, beta_f, alpha_f, gamma, c, v, W, delta),
                             harvest_cournot_own = q_c_w(x, sigma_own, alpha_w, beta_f, alpha_f, gamma, c, v, W, delta),
                             steady_cournot_min = growth - harvest_cournot_min, 
                             steady_cournot_max = growth - harvest_cournot_max, 
                             steady_cournot_own = growth - harvest_cournot_own)

# Looking for the steady state : 

results %>% subset(steady_cournot_min < 0.1 & steady_cournot_min > -0.1)
# One steady state at 91,398 (and also 0)
results %>% subset(steady_cournot_max < 0.1 & steady_cournot_max > -0.1)
# One steady state at 91,370 (and also 0)
results %>% subset(steady_cournot_own < 0.1 & steady_cournot_own > -0.1)
# 3 steady states : 1,343 ; 7,297 ; 91,359

### IV. Betrand Competition ######
# Parameters for demand function (we worked with inverse demand function up to now)
# Seems like there are typos in their functions. Should not matter as long as alpha_w=alpha_f and beta_w=beta_f
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))

b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))

# Cannot make sense of how this works, my function always ends up below the monopoly harvest, while the proof shows it should
# be larger : 
# i. Rewrite the functions from scratch and verify if it works
# ii. If not, look at the derivation of q_b in equation (30)
# iii. If still not : look at the proof, but numerical simulations were consistent with their argument on v_tilde in proof of Lemma
# Z.a


results = results %>% mutate(harvest_bertrand_min = q_b_w(x, b_w, b_f, a_w, a_f, e, sigma_min, v, c, Omega_min),
                             harvest_bertrand_max = q_b_w(x, b_w, b_f, a_w, a_f, e, sigma_max, v, c, Omega_max),
                             harvest_bertrand_own = q_b_w(x, b_w, b_f, a_w, a_f, e, sigma_own, v, c, Omega_own),
                             steady_bertrand_min = growth - harvest_bertrand_min, 
                             steady_bertrand_max = growth - harvest_bertrand_max, 
                             steady_bertrand_own = growth - harvest_bertrand_own)



results %>% subset(steady_bertrand_min > -0.1 & steady_bertrand_min < 0.1)

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour ='growth'))+
  geom_point(aes(y=harvest_bertrand_min, colour = 'harvest_betrand_min'))+
  geom_point(aes(y=harvest_cournot_min, colour = 'harvest_cournot_min'))+
  geom_point(aes(y=harvest_monop_min, colour = 'harvest_monop_min'))

results %>% subset(steady_bertrand_max > -0.1 & steady_bertrand_max < 0.1)
results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour ='growth'))+
  geom_point(aes(y=harvest_bertrand_max, colour = 'harvest_betrand_max'))+
  geom_point(aes(y=harvest_cournot_max, colour = 'harvest_cournot_max'))+
  geom_point(aes(y=harvest_monop_max, colour = 'harvest_monop_max'))

results %>% subset(steady_bertrand_own > -0.1 & steady_bertrand_own < 0.1)
results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour ='growth'))+
  geom_point(aes(y=harvest_bertrand_own, colour = 'harvest_betrand_own'))+
  geom_point(aes(y=harvest_cournot_own, colour = 'harvest_cournot_own'))+
  geom_point(aes(y=harvest_monop_own, colour = 'harvest_monop_own'))

#Could it please work?