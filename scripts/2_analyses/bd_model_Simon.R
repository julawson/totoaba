library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(ggplot2)
library(rootSolve)
#########################################################################################################
# 1. Reprendre le papier avec TOUTES LES PREUVES : on a pas les mêmes résultats qu'eux. 
# 2. Reprendre toutes les justifications: ça m'énerve là putain 
# 3. Il faut simplement reprendre en gardant en tête qu'ils ont beaucoup d'erreurs de calcul et certaines
# de raisonnement dans les preuves. 
########################################################################################################

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

x = seq(0,k,1) # Generate space for x

### II. Monopoly set-up #####

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
harvest_monop = function(x, sigma, alpha, c, beta, W){
  y = ((alpha - c)*(sigma^2)*(x^2))/(2*beta*(sigma^2)*(x^2) + 2*W)
  return(y)
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

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop_own, colour = 'Monop'))+
  geom_point(aes(y=harvest_cournot_own, colour = 'Cournot'))

### IV. Bertrand Competition ######
# Parameters for demand function (we worked with inverse demand function up to now)
# Redid the math : coefficients are good. 
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
# --> Their equation is wrong!!!

# Closed form solution : problem
q_w_b = function(a_f, a_w, b_f, b_w, e, sigma, c, v, W, x){
  num = (sigma^2)*(x^2)*b_w*(b_f * (2*a_w + e*v - 2*b_w*c) + e*(c + a_f))
  den = (sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2))
  y   = num/den
  return(y)
}
# Use price equations and wage to recover quantities
# Price equations from the manuscript are correct
p_f_b = function(a_f, a_w, b_f, b_w, e, v, s, c){
  y = (2*b_w * (a_f + v*b_f) + e*(a_w + (s + c)*b_w))/(4*b_w*b_f - (e^2))
  return(y)
}

p_w_b = function(a_f, a_w, b_f, b_w, e, v, s, c){
  y = (2*b_f * (a_w + (s + c)*b_w) + e*(a_f +v*b_f))/(4*b_w*b_f - (e^2))
  return(y)
}

# I think their equilibrium wage is not however:

s_b = function(a_f, a_w, b_f, b_w, e, v, c, W, sigma, x){
  y = (W*b_w*(b_f*(2*a_w + e*v) + e*a_f + c*((e^2) + 2*b_w*b_f)))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

s_b_own = function(a_f, a_w, b_f, b_w, e, v, c, W, sigma, x){
  y = (2*W*b_w*(b_f*(2*a_w + v*e) + a_f*e + c*((e^2) - 2*b_f*b_w)))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + 2*b_w*W*(2*b_f*b_w - (e^2)))
  return(y)
}

check_equilibrium_wages = data.frame(x) %>% mutate(wage_paper = s_b(a_f, a_w, b_f, b_w, e, v, c, W, sigma_min, x),
                                                  wage_recalc = s_b_own(a_f, a_w, b_f, b_w, e, v, c, W, sigma_min, x),
                                                  wage_diff = ((wage_paper-wage_recalc)/wage_recalc)*100)

check_equilibrium_wages %>% ggplot(aes(x=x))+
  geom_point(aes(y=wage_paper, colour='wage_paper'))+
  geom_point(aes(y=wage_recalc, colour = 'wage_recalc'))

check_equilibrium_wages %>% ggplot(aes(x=x, y= wage_diff)) + geom_point()

# After recalculating the results, their function is wrong, the proof is wrong. 



# Lemma 2: verify if it holds. 
# Focus on MONOPOLY functions : 
# From equation 14 in the paper : 
harvest_monop = function(x, sigma, alpha, c, beta, W){
  y = ((alpha - c)*(sigma^2)*(x^2))/(2*beta*(sigma^2)*(x^2) + 2*W)
  return(y)
}
# Now in Lemma 2a proof p. 11, they use the definition of linked markets, and push it to the extreme where there
# is no longer any substitutability i.e, gamma=0. 
harvest_monop_lemma = function(x, sigma, a_w, b_w, c, W){
  y = ((sigma^2)*(x^2)*(a_w - b_w*c))/(2*(sigma^2)*(x^2) + W*b_w)
  return(y)
}

# In this case : a_m = alpha/beta, b_m = 1/beta
a_m = alpha/beta
b_m = 1/beta

comp_monop = data.frame(x) %>% mutate(growth = growth(x, r, k),
                                      harvest_monop_true = harvest_monop(x, sigma_min, alpha, c, beta, W),
                                      harvest_monop_arg_lemma = harvest_monop_lemma(x, sigma_min, a_m, b_m, c, W), 
                                      harvest_monop_false = harvest_monop_lemma(x, sigma_min, a_w, b_w, c, W),
                                      harvest_monop_true_own = harvest_monop(x, sigma_own, alpha, c, beta, W),
                                      harvest_monop_arg_lemma_own = harvest_monop_lemma(x, sigma_own, a_m, b_m, c, W))

comp_monop %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop_true, colour='True monopoly function'))+
  geom_point(aes(y=harvest_monop_false, colour = 'Wrong function with wrong coeff'))+
  geom_point(aes(y=harvest_monop_arg_lemma, colour = "Wrong monopoly function"))

comp_monop %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop_true_own, colour='True monopoly function'))+
  geom_point(aes(y=harvest_monop_arg_lemma_own, colour = "Wrong monopoly function"))

# Shows that they use the wrong function to do their proof : it it unclear where the Bertrand function sits in that case. 
# Lemma 2 rests on the assumption that a_w = a_m. 
# However, if one wants to use the price derivation for the monopoly, it is important to consider that :
# 1) gamma= 0
# 2) This calls for a new calculation of a_w and a_m, because they are not the same. 


# With renewed computation : 

harvest_bertrand_own = function(a_f, a_w, b_f, b_w, x, sigma, e, v, W, c){
  y = ((sigma^2) * (x^2) * b_w * (b_f * (2*a_w + e*v) + c*((e^2) - 2*b_f*b_w) + e*a_f))/(2*W*b_w * (2*b_f*b_w - (e^2)) + (4*b_w*b_f - (e^2)) * (sigma^2) * (x^2))
  return(y)
}

results = results %>% mutate(harvest_bertrand_min = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma_min, e, v, W, c),
                             harvest_bertrand_max = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma_max, e, v, W, c),
                             harvest_bertrand_own = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma_own, e, v, W, c))
                             
results %>% ggplot(aes(x=x)) +
  geom_point(aes(y=growth, colour ='Growth'))+
  geom_point(aes(y=harvest_monop_own, colour = 'Monopoly'))+
  geom_point(aes(y=harvest_cournot_own, colour = "Cournot"))+
  geom_point(aes(y=harvest_betrand_own, colour = 'Bertrand'))

# We end up with weird results : 
# How could a Bertrand model end up to less harvest than the monopoly? 
# The intuition is to lower the prices. In order to do that, it has to produce more, everything else equal.
# If lower harvest, it means there is a larger price : does not seem possible.

# Need to do a little sensitivity analysis : 
gamma = 1
delta = 4*beta_f*beta_w-gamma^2
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))

b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))



sensitivity = data.frame(x) %>% mutate(growth = growth(x, r, k), 
                                       harvest_monop = harvest_monop(x, sigma_own, alpha, c, beta, W),
                                       harvest_cournot = q_c_w(x, sigma_own, alpha_w, beta_f, alpha_f, gamma, c, v, W, delta), 
                                       harvest_bertrand = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma_own, e, v, W, c))

sensitivity %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop, colour = 'Monop'))+
  geom_point(aes(y=harvest_cournot, colour = 'Cournot'))+
  geom_point(aes(y=harvest_bertrand, colour = 'Bertrand'))
# Weird result : there must still be some error somewhere. With a perfectly substitutable good, should have a larger production. 

gamma = 1
delta = 4*beta_f*beta_w-gamma^2
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))

v = 100
# Set lower costs for acquaculture : does not seem to change the ordering. 

sensitivity = data.frame(x) %>% mutate(growth = growth(x, r, k), 
                                       harvest_monop = harvest_monop(x, sigma_own, alpha, c, beta, W),
                                       harvest_cournot = q_c_w(x, sigma_own, alpha_w, beta_f, alpha_f, gamma, c, v, W, delta), 
                                       harvest_bertrand = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma_own, e, v, W, c))

sensitivity %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop, colour = 'Monop'))+
  geom_point(aes(y=harvest_cournot, colour = 'Cournot'))+
  geom_point(aes(y=harvest_bertrand, colour = 'Bertrand'))

# Reprendre ça : ça ne fonctionne pas comme prévu. 
# 1. Les calculs semblent faux de la part de l'article. 
# 2. Ceci dit, même avec l'argument de l'article ça ne semble pas prendre : 
# Faisons:
gamma = 0.75
a_m = a_w
b_m = b_w

article_res = data.frame(x) %>% mutate(growth = growth(x, r, k), 
                                       monop = harvest_monop_lemma(x, sigma_min, a_w, b_w, c, W), 
                                       cournot = q_c_w(x, sigma_min,alpha_w, beta_f, alpha_f, gamma, c, v, W, delta),
                                       bertrand = q_w_b(a_f, a_w, b_f, b_w, e, sigma_min, c, v, W, x))
article_res %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=monop, colour = 'Monop'))+
  geom_point(aes(y=cournot, colour = 'Cournot'))+
  geom_point(aes(y=bertrand, colour = 'Bertrand'))

# Toujours un problème : la quantité de Cournot est au dessus de celle du monopole. 
# On devrait avoir : Cournot un poil en dessous, Bertrand un poil au dessus. 
# Ma fonction de monopole est juste. 

# Conclusion : 
# Il faut reprendre à partir de Cournot, au moins. 
# Ensuite, il faut reprendre le Bertrand, et les lemmes aussi. Ce n'est vraiment pas clair leur truc. 


######

#############################################################################################################
## Rewrite of the model ######
########################################################
results = data.frame(x)
### I. Monopoly ####
# My computations yield the same results as their
harvest_monop = function(sigma, x, W, alpha, c, beta){
  y = ((sigma^2)*(x^2)*(alpha - c))/(2*(sigma^2)*(x^2)*beta + 2*W)
  return(y)
}
results = results %>% mutate(harvest_monop_min = harvest_monop(sigma_min, x, W, alpha, c, beta),
                             growth = growth(x,r,k),
                             harvest_monop_own = harvest_monop(sigma_own, x, W, alpha, c, beta))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth'))+
  geom_point(aes(y=harvest_monop_own, colour = 'Monop'))

### II. Cournot #####
# Reaction functions : same result
# Closed form depending on s : same result
# Equilibrium poacher wage : different result 

s_c_paper = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (W*(2*beta_f*alpha_w - gamma*(alpha_f - v) - 2*c*beta_f))/((sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)) + 2*W*beta_f)
  return(y)
}

s_c_own = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (2*W*(beta_f*(alpha_w - c) - gamma*(alpha_f + v)))/((sigma^2)*(x^2)*(4*beta_f*beta_w - gamma^2) + 4*beta_f*W)
  return(y)
}

results  = results %>% mutate(wage_cournot_paper = s_c_paper(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x),
                   wage_cournot_own = s_c_own(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=wage_cournot_paper, colour = 'Paper wage'))+
  geom_point(aes(y=wage_cournot_own, colour = 'Recomputed wage'))
# End up with a big magnitude difference, from an over estimation by 160% to 40%:
results %>% ggplot(aes(x=x))+
  geom_point(aes(y=(wage_cournot_own-wage_cournot_paper)/wage_cournot_own, colour = 'Paper wage'))

# Substituting into demand function : they make a mistake, by a factor of 2. 
# Compare both functions : 
q_cournot_wild_paper = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*alpha_w * beta_f - gamma*(alpha_f - v) - 2*beta_f*c))/(2*beta_f*W + (sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)))
  return(y)
}
q_cournot_wild_own = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_f*(alpha_w - c) - gamma*(alpha_f + v)))/(4*beta_f*W + (4*beta_f*beta_w - (gamma^2))*(sigma^2)*(x^2))
  return(y)
}
results = results %>% mutate(q_cournot_wild_og = q_cournot_wild_paper(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own = q_cournot_wild_own(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_og_min = q_cournot_wild_paper(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own_min = q_cournot_wild_own(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=q_cournot_wild_og_min, colour = 'Cournot Paper'))+
  geom_point(aes(y=q_cournot_wild_own_min, colour = 'Cournot own'))+
  geom_point(aes(y=harvest_monop_min, colour = 'Monop'))

sum((results$harvest_monop_own - results$q_cournot_wild_og)>0)
sum((results$harvest_monop_own - results$q_cournot_wild_own)>0)

# Problem with their formulation : they cannot say unambiguously that cournot < monopoly. 
# I propose we move with my formulation : q_cournot_wild_own. 

### III. Bertrand #####
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))
# Not the same result for salary
s_b_paper = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = (W*b_w*(b_f*(2*a_w + e*v) + e*a_f + c*((e^2) + 2*b_w*b_f)))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

s_b_own = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = (2*W*b_w*( b_f*(2*a_w + e*v) + c*(e^2 - 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + 2*W*b_w*(2*b_w*b_f - (e^2)))
  return(y)
}

results = results %>% mutate(s_b_own = s_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x),
                             s_b_paper = s_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=s_b_own, colour = 'Own wage'))+
  geom_point(aes(y=s_b_paper, colour = 'Their wage'))

q_b_own = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v) + c*((e^2)- 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + 2*W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

q_b_paper = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v - 2*b_w*c) + e*(c+a_f)))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

results = results %>% mutate(q_b_paper = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x),
                             q_b_own = q_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x))

### IV. Wrap up their results #####

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=growth, colour = 'Growth')) +
  geom_point(aes(y=harvest_monop_min, colour= 'Monop')) +
  geom_point(aes(y=q_cournot_wild_og, colour = 'Cournot')) +
  geom_point(aes(y=q_b_paper, colour='Bertrand'))

### V. Conclusion #####
# Wrong assumption for their demonstration: cannot assume that a_m  = a_w
# If one does, it's the same as demonstrating that the monopoly produces less than 
# a Bertrand when there is no competition. 
# Or assuming the behavior of a monopoly on a duopolistic market. It is not really comparing both outcomes

# Their theoretical proofs are wrong : using their values and functions, we do not find q_m > q_b
results %>% ggplot(aes(x=x))+
  geom_line(aes(y = growth, colour = 'growth'))+
  geom_line(aes(y = harvest_monop_min, colour = 'Monop'))+
  geom_line(aes(y = q_cournot_wild_og_min, colour ='Cournot'))+
  geom_line(aes(y = q_b_paper, colour = 'Bertrand'))
# This result is difficult to interpret : 
# 1. Cournot does result in a smaller harvest, and can have better results
# 2. Bertrand would allegedly lead to a large steady state as well. However, with those results
# it seems like the marginal cost curve of the poachers is too high : s+c is huge, way larger then 
# v. Therefore, the only steady state to arise would be 0. 
# 3. It would be interesting to have a fresh set of eyes on this : I think their result lies on the fact that
# they assume *a_m= a_w* to prove what they do, and they shouldn't. Moreover, their results even with their formulations
# do depend on a non-negative v (do I need to prove they're wrong mathematically? The simulations show they're wrong.)
# 4. Indeed, there's a range of values where it makes sense to say that bertrand leads to more harvesting : when it's
# more profitable for them to do so. But there exists a point, considering the nature of the resource, 
# such that it costs too much money. 

### VI. Sensitivity analysis#####
# 
q_b_500 = ((sigma_own^2)*(x^2)*b_w*(b_f*(2*a_w + e*500 - 2*b_w*c) + e*(c + a_f)))/((sigma_own^2)*(x^2)*(4*b_w*b_f - (e^2)) +W*b_w*(2*b_f*b_w - (e^2)))
q_b_1000 = ((sigma_own^2)*(x^2)*b_w*(b_f*(2*a_w + e*v - 2*b_w*c) + e*(c + a_f)))/((sigma_own^2)*(x^2)*(4*b_w*b_f - (e^2)) +W*b_w*(2*b_f*b_w - (e^2)))
q_b_5000 = ((sigma_own^2)*(x^2)*b_w*(b_f*(2*a_w + e*5000 - 2*b_w*c) + e*(c + a_f)))/((sigma_own^2)*(x^2)*(4*b_w*b_f - (e^2)) +W*b_w*(2*b_f*b_w - (e^2)))
q_b_10000 = ((sigma_own^2)*(x^2)*b_w*(b_f*(2*a_w + e*10000 - 2*b_w*c) + e*(c + a_f)))/((sigma_own^2)*(x^2)*(4*b_w*b_f - (e^2)) +W*b_w*(2*b_f*b_w - (e^2)))

monop = ((alpha -c)*(sigma_own^2)*(x^2))/(2*beta_w *(sigma_own^2)*(x^2)+2*W)
ultime = data.frame(x, q_b_1000, q_b_5000, q_b_10000,q_b_500, monop, growth = growth(x, r, k))
ultime %>% ggplot(aes(x=x))+
  geom_line(aes(y=q_b_1000, colour="q_b : v=1000"), size = 1)+
  geom_line(aes(y=q_b_5000, colour="q_b : v=5000"), size = 1)+
  geom_line(aes(y=q_b_500, colour="q_b : v=500"), size = 1)+
  geom_line(aes(y=monop, colour = "Monop"), size = 1)+
  geom_line(aes(y=growth, colour = 'Growth'), size = 1)
# This shows that the results are indeed depending on v. 

ultime = ultime %>% mutate(q_b_own_500 = q_b_own(a_f, a_w, b_f, b_w, c, e, 500, W, sigma_own, x), 
                           q_b_own_2000 = q_b_own(a_f, a_w, b_f, b_w, c, e, 2000, W, sigma_own, x), 
                           q_b_own_5000 = q_b_own(a_f, a_w, b_f, b_w, c, e, 5000, W, sigma_own, x), 
                           q_b_own_1000 = q_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x))
ultime %>% ggplot(aes(x=x))+
  #geom_line(aes(y=q_b_own_500, colour="q_b : v=500"), size = 1)+
  #geom_line(aes(y=q_b_own_2000, colour="q_b : v=2000"), size = 1)+
  geom_line(aes(y=q_b_own_1000, colour="q_b : v=1000"), size = 1)+
  #geom_line(aes(y=q_b_own_5000, colour="q_b : v=5000"), size = 1)+
  geom_line(aes(y=monop, colour = "Monop"), size = 1)+
  geom_line(aes(y=growth, colour = 'Growth'), size = 1)

# My approach to lemma 2 : 
# 
# Condition toujours : 

rhs = ((a_m - b_m *c)/(b_m*(2 * sigma^2 * x^2 + 2*W*b_m)))*(sigma^2 * x^2 * (4*b_f*b_w - e^2) + 2*W*b_w*(2*b_f*b_w - e^2))
rhs = (rhs - e*a_f + c*(2*b_f*b_w - e^2) - 2*b_f*a_w)/(2*b_f*e)
sp = data.frame(x, rhs, f) %>% ggplot(aes(x=x, y=rhs))+geom_line()
sp+geom_hline(yintercept = v)+geom_hline(yintercept = 4000)

ultime%>%ggplot(aes(x=x))+
  geom_line(aes(y=monop, colour='Monop'))+
  geom_line(aes(y=q_b_own(a_w, a_f, b_w, b_f, W, e, v, c, x, sigma), colour='Bertrand v=1000'))+
  geom_line(aes(y=q_b_own(a_w, a_f, b_w, b_f, W, e, 6000, c, x, sigma), colour='Bertrand v=4000'))+
  geom_vline(xintercept = 11000)

# Further along the proof of lemma 2: 
# Redefine their function : 
q_w_b_paper = function(sigma, x, a_f, a_w, b_f, b_w, c, v, e, W){
  y = (sigma^2 * x^2 * b_w * (b_f*(2*a_w + e*v -2*b_w*c) + e*(c+a_f)))/(sigma^2 * x^2 * (4*b_f * b_w - e^2) + W*b_w*(2*b_f*b_w - e^2))
  return(y)
}

q_monop = function(a_m, b_m, c, sigma, x, W){
  y = sigma^2 * x^2*(a_m - b_m*c)/(2*(sigma^2)*(x^2) +2* W*b_m)
  return(y)
}

harvest_monop_q = function(alpha, beta, c, W, sigma, x){
  y = (sigma^2 * x^2 * (alpha - c))/(2 * sigma^2 * x^2 * beta + 2*W)
}

data.frame(x, 
           monop_aw = q_monop(a_w, b_w, c, sigma, x, W), 
           bertrand = q_w_b_paper(sigma, x, a_f, a_w, b_f, b_w, c, v, e, W),
           bertrand_5000 = q_w_b_paper(sigma, x, a_f, a_w, b_f, b_w, c, 5000, e, W),
           monop_am = q_monop(a_m, b_m, c, sigma, x, W), 
           monop_true = harvest_monop_q(alpha, beta, c, W, sigma, x)) %>%
  ggplot(aes(x))+
  geom_line(aes(y=monop_aw, colour = 'Monop with a_w'), size=1.5)+
  geom_line(aes(y=bertrand, colour= 'Bertrand, v=1000'), size=1.5)+
  geom_line(aes(y=bertrand_5000, colour= 'Bertrand, v=5000'), size=1.5)+
  geom_line(aes(y=monop_am, colour = 'Monop with a_m'), size=1.5)+
  geom_line(aes(y=monop_true, colour = 'True monopoly function'), size=1.5)
# THis is the thing to remember : 
# They proved a trivial result : a monopoly produces more than a Bertrand when there is no competition.
# However, it becomes more complicated as things move forward : the results are not cost independent. 
# One cannot restrict the analysis to the case where it is guaranteed that there are two competitors. 

# One question is how do they get to a large equilibrium in theory, like how would it work? 

# How do my calculations fare in that case? 
gamma = 0.75
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))
data.frame(x, 
           monop_aw = q_monop(a_w, b_w, c, sigma, x, W), 
           bertrand = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, v, W, c),
           bertrand_5000 = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, 5000, W, c),
           bertrand_3400 = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, 3400, W, c),
           
           monop_am = q_monop(a_m, b_m, c, sigma, x, W), 
           monop_true = harvest_monop_q(alpha, beta, c, W, sigma, x)) %>%
  ggplot(aes(x))+
  geom_line(aes(y=monop_aw, colour = 'Monop with a_w'), size=1.5)+
  geom_line(aes(y=bertrand, colour= 'Bertrand, v=1000'), size=1.5)+
  geom_line(aes(y=bertrand_5000, colour= 'Bertrand, v=5000'), size=1.5)+
  geom_line(aes(y=bertrand_3400, colour= 'Bertrand, v=3400'), size=1.5)+
  geom_line(aes(y=monop_am, colour = 'Monop with a_m'), size=1.5)+
  geom_line(aes(y=monop_true, colour = 'True monopoly function'), size=1.5)
# Interesting : if we believe my calculations are better, there is no more crossing. 
# In this case, depending on the value of the marginal cost of producing for the aquaculture sector, 
# there is gonna be a lower equilibrium, or an unlikely larger equilibrium. I'd need to talk about this with someone
# in terms of economics. 
# Crossing is driven by how they systematically seem to forget a 2 on the denominator.
gamma = 0.2
alpha_f = 8000
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_w - (gamma^2))

data.frame(x, 
           monop_aw = q_monop(a_w, b_w, c, sigma, x, W), 
           bertrand = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, v, W, c),
           bertrand_5000 = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, 7000, W, c),
           bertrand_3400 = harvest_bertrand_own(a_f, a_w, b_f, b_w, x, sigma, e, 3400, W, c),
           
           monop_am = q_monop(a_m, b_m, c, sigma, x, W), 
           monop_true = harvest_monop_q(alpha, beta, c, W, sigma, x)) %>%
  ggplot(aes(x))+
  geom_line(aes(y=monop_aw, colour = 'Monop with a_w'), size=1.5)+
  geom_line(aes(y=bertrand, colour= 'Bertrand, v=1000'), size=1.5)+
  geom_line(aes(y=bertrand_5000, colour= 'Bertrand, v=5000'), size=1.5)+
  #geom_line(aes(y=bertrand_3400, colour= 'Bertrand, v=3400'), size=1.5)+
  geom_line(aes(y=monop_am, colour = 'Monop with a_m'), size=1.5)+
  geom_line(aes(y=monop_true, colour = 'True monopoly function'), size=1.5)

# The degree of substitutability is going to be key in this result : it does change the harvest results for the Bertrand equilibrium
# It pushed the harvest rates up : makes sense. 
# May need to make it endogenous : it entails significant differences in results. 

# How to make profit function depend on choice of gamma, and remain linear in quantities and prices. 
