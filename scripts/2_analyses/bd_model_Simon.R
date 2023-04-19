library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(ggplot2)
library(rootSolve)
####################################################################################################################
# Model for Damania, Bulte - The economics of wildlife farming and endangered species conservation - 2005
####################################################################################################################

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

### II. Perfect competition ####
# Define growth function

growth <- function(x,r,k){
  return(r*x*(1-x/k))
}

#Initiate results dataframe
results = data.frame(x)

results$growth = growth(x,r,k)  

harvest_competition = function(x, sigma){
  y = P*sigma^2*x^2/(4*W)
  return(y)
}


### III. Monopoly set-up #####

# Define harvest function in Monopoly : 
# - Computations yield the same results.
harvest_monop = function(x, sigma, alpha, c, beta, W){
  y = ((alpha - c)*(sigma^2)*(x^2))/(2*beta*(sigma^2)*(x^2) + 2*W)
  return(y)
}
# Set results 
results$harvest_monop_min = harvest_monop(x, sigma_min, alpha, c, beta, W)
results$harvest_monop_max = harvest_monop(x, sigma_max, alpha, c, beta, W)
results$harvest_monop_own = harvest_monop(x, sigma_own, alpha, c, beta, W)

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
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=growth, colour = 'Growth'))+
  geom_line(aes(y=harvest_monop_min, colour = 'Monop. min'))+
  geom_line(aes(y=harvest_monop_own, colour = 'Monop. own'))

### IV. Cournot #####
# Parameters: 
    # Demand
beta_f = beta
beta_w = beta
alpha_w = alpha
alpha_f = alpha
gamma = 0.75
    # Costs:
v = 1000
# Reaction functions : same result
# Closed form depending on s : same result
# Equilibrium poacher wage : different result 

s_c_paper = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (W*(2*beta_f*alpha_w - gamma*(alpha_f - v) - 2*c*beta_f))/((sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)) + 2*W*beta_f)
  return(y)
}

s_c_own = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (2*W*(2*beta_f*(alpha_w - c) - gamma*(alpha_f + v)))/((sigma^2)*(x^2)*(4*beta_f*beta_w - gamma^2) + 4*beta_f*W)
  return(y)
}

results  = results %>% mutate(wage_cournot_paper = s_c_paper(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x),
                   wage_cournot_own = s_c_own(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=wage_cournot_paper, colour = 'Paper wage'))+
  geom_point(aes(y=wage_cournot_own, colour = 'Recomputed wage'))
# End up with a big magnitude difference, and a crossing.


# Substituting into demand function : they make a mistake, by a factor of 2 : say it's a typo
# Compare both functions : 
q_cournot_wild_paper = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*alpha_w * beta_f - gamma*(alpha_f - v) - 2*beta_f*c))/(2*beta_f*W + (sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)))
  return(y)
}
q_cournot_wild_own = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/(4*beta_f*W + (4*beta_f*beta_w - (gamma^2))*(sigma^2)*(x^2))
  return(y)
}
results = results %>% mutate(q_cournot_wild_og = q_cournot_wild_paper(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own = q_cournot_wild_own(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_og_min = q_cournot_wild_paper(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own_min = q_cournot_wild_own(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W))

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=q_cournot_wild_og, colour = 'Cournot Paper'))+
  geom_point(aes(y=q_cournot_wild_own, colour = 'Cournot own'))
# Their harvest overestimates actual harvest

results %>% ggplot(aes(x=x))+
  geom_point(aes(y=q_cournot_wild_og, colour = 'Cournot Paper'))+
  geom_point(aes(y=harvest_monop_own, colour = 'Monop'))+
  geom_point(aes(y= q_cournot_wild_own, colour = 'Cournot recalc'))
# From this graph : 
# Their function does not yield an unambiguously lower harvest for Cournot, our does:
sum((results$harvest_monop_own-results$q_cournot_wild_og)>0)/k
sum((results$harvest_monop_own-results$q_cournot_wild_own)>0)/k


results %>% mutate(ss_cournot_own = growth - q_cournot_wild_own) %>%
  subset(ss_cournot_own > -0.1 & ss_cournot_own < 0.1)
# With our result, we only have a high steady state at 91370, with sigma = sigma_own

results %>% mutate(ss_cournot_own = growth - q_cournot_wild_og) %>%
  subset(ss_cournot_own > -0.1 & ss_cournot_own < 0.1)
# With their result : 3 steady states : 
# - 1344
# - 7298
# - 91360

# Recompute the farmer's production function : 
q_farmed_cournot = function(alpha_w, alpha_f, beta_w, beta_f, s, c, v, gamma){
  y = (2*beta_w * (alpha_f - v) - gamma*(alpha_w - s - c))/(4*beta_w*beta_f - gamma^2)
}
results  = results %>% mutate(q_farmed_cournot(alpha_w, alpha_f, beta_w, beta_f, s_c_own(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_own, x), c, v, gamma))

# Conclusion : 
# They have a typo in their result, and do not use the reported sigma in the documents they cite. 
# However, with q_cournot_own, and sigma_own we find the same results

### V. Bertrand #####
# Define parameters for inverse demand function : 
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_f - (gamma^2))

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
# Our wage is larger, and crosses, again. 

q_b_own = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v) + c*((e^2)- 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + 2*W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

q_b_paper = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v - 2*b_w*c) + e*(c+a_f)))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

results = results %>% mutate(q_b_paper = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x),
                             q_b_paper_min = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x),
                             q_b_own = q_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x))
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=q_b_own, colour = 'Closed form own'))+
  geom_line(aes(y=q_b_paper, colour = 'Paper'))

# Find the steady states : 
# Reference in the paper says 0; 8,480 ; 83,470

results %>% 
  mutate(ss_bertrand = growth - q_b_paper)%>%
  subset(ss_bertrand > -.1 & ss_bertrand <.1)


results %>% 
  mutate(ss_bertrand = growth - q_b_paper_min)%>%
  subset(ss_bertrand > -.1 & ss_bertrand <.1)

results %>% 
  mutate(ss_bertrand = growth - q_b_own)%>%
  subset(ss_bertrand > -.1 & ss_bertrand <.1)

# I cannot find ways to replicate their findings.


### VI. Wrap up the results #####

# This graph summarizes their findings following the paper's function
results %>% ggplot(aes(x = x))+
  geom_point(aes(y = growth, colour = 'Growth')) +
  geom_point(aes(y = harvest_monop_own, colour= 'Monop')) +
  geom_point(aes(y = q_cournot_wild_og, colour = 'Cournot')) +
  geom_point(aes(y = q_b_paper, colour='Bertrand'))

# With our results : 
results %>% ggplot(aes(x = x)) +
  geom_point(aes(y = growth, colour = 'Growth'))+
  geom_point(aes(y = harvest_monop_own, colour = 'Monop'))+
  geom_point(aes(y = q_cournot_wild_own, colour = 'Cournot'))+
  geom_point(aes(y = q_b_own, colour = 'Bertrand'))

# I am confused because : 
# 1. I cannot seem to rationalize their findings for the Bertrand steady state.
# 2. With their function, Lemma 2 does not stand: the poaching level in the case of Bertrand competition
# is not larger than in the case of a monopoly

### VII. Rationalizing the way they work for Lemma 2 #####

# Wrong assumption for their demonstration: 

# They assume that : a_m  = a_w
# These coefficients are computed using the assumption of two producers on the market. 
# Therefore, using a_w would mean : 
# - There is indeed a competitor on the market, and the firm still behaves like a monopoly. 
# In that case, it would mean that the firm behaves like a monopolist on a duopolistic market. Because it 
# operates with residual demand, it produces less. 

# Check : https://www.overleaf.com/5228619579rgmypmgnvbmn

# Therefore, the proof does not compare the right counterfactuals. Let's illustrate that : 

# A. Illustrate the difference between a monopoly operating on a monopolistic and duopolistic market
harvest_monop_indirect = function(a_m, b_m, c, W, sigma, x){
  y = (a_m - b_m *c)*(sigma^2)*(x^2)/(2*sigma^2*x^2 + 2*b_m*W)
  return(y)
}
a_m = alpha/beta
b_m = 1/beta


lemma = data.frame(x, 
                   growth = growth(x,r,k),
                   monop_aw = harvest_monop_indirect(a_w, b_w, c, W, sigma_min, x),
                   monop_am = harvest_monop_indirect(a_m, b_m, c, W, sigma_min, x),
                   q_b_aw = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x))

colors = c('blue', 'grey', 'black')
lemma %>% ggplot(aes(x=x))+
  geom_line(aes(y = monop_aw, colour = 'Monopoly on a duopolistic market'))+
  geom_line(aes(y = monop_am, colour = 'Monopoly on a monopolistic market'))+
  geom_line(aes(y = q_b_aw, colour = 'Bertrand on a duopolistic market'))+
  scale_color_manual(values=colors)+
  xlab('x')+
  ylab("Harvest level")+
  ggtitle('Comparison of Bertrand and Monopoly harvests depending on market structure')+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="bottom")
  
# This result shows that if we do not (wrongfully so) assume that a_m = a_w, then their result
# do not hold anymore. 
colors = c('blue', 'grey', 'black', 'red')

# Out of curiosity:  
lemma %>%
  mutate(q_cournot_wild = q_cournot_wild_paper(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W))%>%
  ggplot(aes(x=x))+
  geom_line(aes(y = monop_aw, colour = 'Monopoly with a_w'))+
  geom_line(aes(y = monop_am, colour = 'Monopoly with a_m'))+
  geom_line(aes(y = q_b_aw, colour = 'Bertrand'))+
  geom_line(aes(y = q_cournot_wild, colour = 'Cournot'))
  scale_color_manual(values=colors)
# This result shows the inconsistency in their method and results. 

colors = c('blue', 'grey', 'black', 'red')
lemma = data.frame(x=seq(0,12500), 
                   growth = growth(seq(0,12500),r,k),
                   monop_aw = harvest_monop_indirect(a_w, b_w, c, W, sigma_own, seq(0,12500)),
                   monop_am = harvest_monop_indirect(a_m, b_m, c, W, sigma_own, seq(0,12500)),
                   q_b_aw = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, seq(0,12500)))
# Out of curiosity:  
lemma %>%
  mutate(q_cournot_wild = q_cournot_wild_paper(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W))%>%
  ggplot(aes(x=x))+
  geom_line(aes(y = monop_aw, colour = 'Monopoly with a_w'))+
  geom_line(aes(y = monop_am, colour = 'Monopoly with a_m'))+
  geom_line(aes(y = q_b_aw, colour = 'Bertrand'))+
  geom_line(aes(y = growth, colour = 'Growth'))+
  scale_color_manual(values=colors)+
  xlab('x')+
  ylab("Harvest level")+
  ggtitle('Comparison of Bertrand and Monopoly steady states on market structure')+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="bottom")
