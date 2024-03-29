rm(list=ls())
library(here)
library(tidyverse)
library(janitor)
library(broom)
library(deSolve)
library(ggplot2)
library(rootSolve)

####################################################################################################################
# Model for Damania, Bulte - The economics of wildlife farming and endangered species conservation - 2007
####################################################################################################################

### I. Parameters #####
# Define the parameters for the biological resource
r <- 0.16 #intrinsic rate of population increase, estimated
k <- 100000 #carrying capacity, number of individuals
sigma_min <- 0.000167 #catchability, ranging from 0.000167 to 0.000256
sigma_max <- 0.000256
sigma_own = 0.000367718

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

### II. Perfect competition : the basic poaching model ####

# Define growth function
growth <- function(x,r,k){
  return(r*x*(1-x/k))
}

#Initiate results dataframe
results = data.frame(x)%>%
  mutate(growth = growth(x,r,k))

# Define poaching in this scenario
harvest_competition = function(x, sigma){
  y = P*sigma^2*x^2/(4*W)
  return(y)
}


### III. Monopoly set-up #####

# Define harvest function in Monopoly : 
# - Mathematical computations yield the same results.
harvest_monop = function(x, sigma, alpha, c, beta, W){
  y = ((alpha - c)*(sigma^2)*(x^2))/(2*beta*(sigma^2)*(x^2) + 2*W)
  return(y)
}
# Set results 
# Compute growth - harvest
results = results %>% mutate(harvest_monop_min = harvest_monop(x, sigma_min, alpha, c, beta, W),
                             harvest_monop_max = harvest_monop(x, sigma_max, alpha, c, beta, W),
                             harvest_monop_own = harvest_monop(x, sigma_own, alpha, c, beta, W),
                             steady_monop_min = growth - harvest_monop_min,
                             steady_monop_max = growth - harvest_monop_max,
                             steady_monop_own = growth - harvest_monop_own)

# Find where growth = harvest to find steady states (sign change is necessary)
results %>% 
  subset(steady_monop_min>-0.1 & steady_monop_min<0.1) %>%
  select(x, steady_monop_min)
# Steady state is at 90,120 in this case

results %>%
  subset(steady_monop_max>-0.1 & steady_monop_max<0.1) %>%
  select(x, steady_monop_max)

# Steady state is at 90,054

results %>% 
  subset(steady_monop_own>-0.1 & steady_monop_own<0.1) %>%
  select(x, steady_monop_own)
# 3 steady states : 
# 2,625; 
# 7,346; 
# 90,029

# Values in the paper are : 
# 2,625
# 7,346
# 90,030

# Plot the different results for different sigma values
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=growth, colour = 'Growth'), linewidth = 1.1)+
  geom_line(aes(y=harvest_monop_min, colour = 'Monop. min'), linewidth = 1.1)+
  geom_line(aes(y=harvest_monop_own, colour = 'Monop. own'), linewidth = 1.1)+
  theme_bw()+
  ggtitle('Comparing Monopoly harvest for different sigma values')

### IV. Imperfect competition : Cournot duopoly #####
###### A. Parameters: #####
    # Demand
beta_f = beta
beta_w = beta
alpha_w = alpha
alpha_f = alpha
gamma = 0.75
    # Costs:
v = 1000
###### B. Functions derived in the model : #####
### Reaction functions : same result
### Closed form depending on s : same result
### Equilibrium poacher wage : different result 

# Equilibrium poacher wage in the original model
s_c_paper = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (W*(2*beta_f*alpha_w - gamma*(alpha_f - v) - 2*c*beta_f))/((sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)) + 2*W*beta_f)
  return(y)
}
# Equilibrium poacher wage recomputed
s_c_own = function(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma, x){
  y = (2*W*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/((sigma^2)*(x^2)*(4*beta_f*beta_w - gamma^2) + 4*beta_f*W)
  return(y)
}

s_c_own_use = s_c_own(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_own, x)

# Store results
results  = results %>% mutate(wage_cournot_paper = s_c_paper(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x),
                   wage_cournot_own = s_c_own(alpha_w, alpha_f, beta_f, beta_w, gamma, W, c, v, sigma_min, x))


## Plot the different wages
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=wage_cournot_paper, colour = 'Paper wage'), linewidth=1.2)+
  geom_line(aes(y=wage_cournot_own, colour = 'Recomputed wage'), linewidth=1.2)+
  ylab("Wage")+
  xlab('Stock (x)')+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()
# End up with a big magnitude difference

###### C. Supply and demand functions : ##########
# Not the same mathematical results

# Quantity supplied by traders in the article
q_cournot_wild_paper = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/(2*beta_f*W + (sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)))
  return(y)
}
# Quantity supplied recomputed
q_cournot_wild_own = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_f*(alpha_w - c) - gamma*(alpha_f - v)))/(4*beta_f*W + (sigma^2)*(x^2)*(4*beta_f*beta_w - (gamma^2)))
  return(y)
}

# Quantity farmed recomputed
q_cournot_farmed_own = function(sigma, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W){
  y = ((sigma^2)*(x^2)*(2*beta_w*(alpha_f - v) - gamma*(alpha_w-c)) + 2*W*(alpha_f-v))/((4*beta_f*beta_w - gamma^2)*(sigma^2)*(x^2)+ 4*W*beta_f)
  return(y)
}

# Recompute the farmer's production function with wage: 
q_farmed_cournot = function(alpha_w, alpha_f, beta_w, beta_f, s, c, v, gamma){
  y = (2 * beta_w * (alpha_f - v) - gamma* (alpha_w - (s + c)) )/(4*beta_w*beta_f - gamma^2)
}

# Store results
results  = results %>% mutate(q_farmed_cournot_with_wage = q_farmed_cournot(alpha_w, alpha_f, beta_w, beta_f, s_c_own_use, c, v, gamma))

results = results %>% mutate(q_cournot_wild_og = q_cournot_wild_paper(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own = q_cournot_wild_own(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_og_min = q_cournot_wild_paper(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W),
                             q_cournot_wild_own_min = q_cournot_wild_own(sigma_min, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W), 
                             q_cournot_farmed_own = q_cournot_farmed_own(sigma_own, x, alpha_f, alpha_w, beta_f, beta_w, gamma, c, v, W))

#Graph : differences in quantity functions
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=q_cournot_wild_own, colour = 'Own quantity'),  linewidth=1.3)+
  geom_line(aes(y=q_cournot_wild_og, colour = 'Quantity of paper'),linewidth=1.3)+
  ylab("Quantities")+
  xlab('Stock (x)')+
  scale_color_manual(values = c('blue','red'))+
  theme_bw()
# Their harvest overestimates actual harvest

## Graph : differences in quantity functions and comparison with monopoly harvest
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=q_cournot_wild_og, colour = 'Cournot Paper'), linewidth = 1.1)+
  geom_line(aes(y=harvest_monop_own, colour = 'Monopoly'), linewidth = 1.1)+
  geom_line(aes(y= q_cournot_wild_own, colour = 'Cournot recomputed'), linewidth = 1.1)+
  theme_bw()
# From this graph : 
# Their function does not yield an unambiguously lower harvest for Cournot, our does:
sum((results$harvest_monop_own-results$q_cournot_wild_og)>0)/k #share of article results where monopoly > Cournot
sum((results$harvest_monop_own-results$q_cournot_wild_own)>0)/k #share of own results where monopoly > Cournot


###### D. Steady state results ######
results %>% mutate(ss_cournot_own = growth - q_cournot_wild_own) %>%
  subset(ss_cournot_own > -0.1 & ss_cournot_own < 0.1)%>%
  select(x, ss_cournot_own)
# With our result, we only have a high steady state at 91370, with sigma = sigma_own

results %>% mutate(ss_cournot_own = growth - q_cournot_wild_og) %>%
  subset(ss_cournot_own > -0.1 & ss_cournot_own < 0.1)%>%
  select(x, ss_cournot_own)

# With their result : 3 steady states : 
# - 1344
# - 7298
# - 91360



# Conclusion : 
# They have a typo in their result, and do not use the reported sigma in the documents they cite. 
# However, with q_cournot_own, and sigma_own we find the same results

###### E. Investigation of Lemma 1 ######
######## i. Revisit Lemma 1 of the article #####
# After recomputation of their proof : 
# Define function for the left hand side of the proof of Lemma 1 (section 4.2.1 in doc)
LHS = function(x, sigma){
  z = alpha_f - (1/gamma)*(alpha - c)*((gamma^2 * sigma^2 * x^2 + 2*beta_f*W)/(2*beta_w*sigma^2 * x^2 + 2*W))
  return(z)
}

# Find when v = LHS
data.frame(x, LHS = LHS(x, sigma_own))%>%
  mutate(diff = v - LHS)%>%
  subset(diff < 1 & diff >-1)

# Plot when the condition is no longer met
data.frame(x, LHS = LHS(x, sigma_own))%>%
  ggplot(aes(x=x, y= LHS))+
  geom_line(colour = 'gray', linewidth=1.1)+
  geom_hline(yintercept = v, linetype = "dotted", linewidth=1.1)+
  geom_text(aes( x=1, y=v-150, label = 'y=v'))+
  ylab('LHS')+
  geom_text(aes( x = 7194, y=-5000, label = 'X=7194'))+
  geom_vline(xintercept = 7194, linetype = 'dotted', linewidth = 1.1)+
  theme_bw()

# Lemma 1 does not hold.

# Plot the point where the two curves meet. 
results %>%
  ggplot(aes(x=x))+
  geom_line(aes(y= harvest_monop_own, colour = 'Monopoly'), linewidth = 1.1)+
  geom_line(aes(y = q_cournot_wild_og, colour = 'Cournot'), linewidth = 1.1)+
  geom_vline(xintercept = 7194, linewidth=1.1, linetype = 'dotted')+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  ylab('Quantity')

######## ii. Fixing Lemma 1 ####
# My version of the condition
LHS2 = function(x, sigma){
  y = alpha_f - gamma*(alpha - c)*sigma^2 * x^2 /(2*beta * sigma^2 * x^2 + 2*W)
  return(y)
}

# Find if there is a point where condition is not met
data.frame(x, LHS = LHS2(x, sigma_own))%>%
  mutate(diff = v - LHS)%>%
  subset(diff < 1 & diff >1)

# Plot the condition
data.frame(x, LHS= LHS2(x, sigma_own))%>%
  ggplot(aes(x=x, y= LHS))+
  geom_line(colour = 'gray', linewidth=1.1)+
  geom_hline(yintercept = v, linetype = "dotted", linewidth=1.1)+
  geom_text(aes( x=1, y=v-150, label = 'y=v'))+
  ylab('LHS')+
  geom_vline(xintercept = 7194, linetype = 'dotted', linewidth = 1.1)+
  theme_bw()

# Plot our recomputed function compared with monopoly
results %>%
  ggplot(aes(x=x))+
  geom_line(aes(y= harvest_monop_own, colour = 'Monopoly'), linewidth = 1.1)+
  geom_line(aes(y = q_cournot_wild_own, colour = 'Cournot'), linewidth = 1.1)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  ylab('Quantity')






### V. Imperfect competition : Bertrand duopoly #####
###### A. Parameters for inverse demand function #####
e   = gamma/(beta_w*beta_f - (gamma^2))
a_f = (alpha_f*beta_w - alpha_w*gamma)/(beta_w*beta_f - (gamma^2))
a_w = (alpha_w*beta_f - alpha_f*gamma)/(beta_w*beta_f - (gamma^2))
b_f = beta_f/(beta_w*beta_f - (gamma^2))
b_w = beta_w/(beta_w*beta_f - (gamma^2))

######## B. Recalculate functions and compare with manuscript
######## i. Wage #####

# Wage in manuscript
s_b_paper = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = (W*b_w*(b_f*(2*a_w + e*v) + e*a_f + c*((e^2) + 2*b_w*b_f)))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

# recomputed wage
s_b_own = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = (2*W*b_w*( b_f*(2*a_w + e*v) + c*((e^2) - 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + 2*W*b_w*(2*b_w*b_f - (e^2)))
  return(y)
}

results = results %>% mutate(s_b_own = s_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x),
                             s_b_paper = s_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x))


# Graph of wage differences
results %>% ggplot(aes(x=x))+
  geom_line(aes(y=s_b_own, colour = 'Paper wage'), linewidth=1.1)+
  geom_line(aes(y=s_b_paper, colour = 'Recomputed wage'), linewidth=1.1)+
  theme_bw()+
  ylab('Wage in Bertrand Competition')+
  scale_color_manual(values = c('blue', 'red'))
# Our wage is larger, and crosses, again. 


######## ii. Quantity #####
# Quantity from manuscript
q_b_paper = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v - 2*b_w*c) + e*(c+a_f)))/((sigma^2)*(x^2)*(4*b_f*b_w - (e^2)) + W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}

# Recomputed quantity
q_b_own = function(a_f, a_w, b_f, b_w, c, e, v, W, sigma, x){
  y = ((sigma^2)*(x^2)*b_w*(b_f*(2*a_w+e*v) + c*((e^2)- 2*b_f*b_w) + e*a_f))/((sigma^2)*(x^2)*(4*b_w*b_f - (e^2)) + 2*W*b_w*(2*b_f*b_w - (e^2)))
  return(y)
}


# Verify results with wage and formula as a function of wage
s_b = s_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x)
qty = b_w*(2*b_f*a_w + e*(a_f + v*b_f) + (s_b+c)*(e^2 - 2*b_w*b_f))/(4*b_w*b_f - e^2)

# Quantity farmed : 
quantity_farmed = function(s){
  y = b_f *(2*b_w*a_f + v*(e^2 - 2*b_w*b_f) + e*(a_w + (s+c)*b_w))/(4*b_f*b_w - e^2)
  return(y)
}

results = results %>% mutate(q_b_paper = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x),
                             q_b_paper_min = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x),
                             q_b_own = q_b_own(a_f, a_w, b_f, b_w, c, e, v, W, sigma_own, x),
                             farmed_b = quantity_farmed(s_b))

results %>% 
  mutate(qty_2 = qty,)%>%
  ggplot(aes(x=x))+
  geom_line(aes(y=q_b_own, colour = 'Closed form own'))+
  geom_line(aes(y=q_b_paper, colour = 'Paper'))+
  geom_line(aes(y=qty_2, colour = 'Other method'))+
  ggtitle('Bertrand harvest function')

# Compare monopoly harvest and Betrand harvest
results %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y= harvest_monop_own, colour = 'Monopoly harvest'), linewidth = 1.1)+
  geom_line(aes(y= q_b_paper, colour = 'Bertrand wild harvest'), linewidth = 1.1)+
  theme_bw()+
  scale_color_manual(values = c('green','blue'))+
  ylab('Harvest')

# Plot farmed and poached quantities:
results %>%
  ggplot(aes(x=x))+
  geom_line(aes( y = q_b_own, colour = 'Poached'))+
  geom_line(aes(y = farmed_b, colour = 'Farmed'))+
  geom_line(aes(y = growth, colour = 'Growth'))+
  theme_bw()
###### C.  Find the steady states : ~################
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

# This graph summarizes the findings following the paper's function
results %>% ggplot(aes(x = x))+
  geom_line(aes(y = growth, colour = 'Growth'), linewidth = 1.1) +
  geom_line(aes(y = harvest_monop_own, colour= 'Monop'), linewidth = 1.1) +
  geom_line(aes(y = q_cournot_wild_og, colour = 'Cournot'), linewidth = 1.1) +
  geom_line(aes(y = q_b_paper, colour = 'Bertrand'), linewidth = 1.1)+
  theme_bw()+
  scale_color_manual(values = c('green', 'red', 'orange', 'blue'))

# With our results : 
results %>% ggplot(aes(x = x)) +
  geom_line(aes(y = growth, colour = 'Growth'), linewidth = 1.1)+
  geom_line(aes(y = harvest_monop_own, colour = 'Monop'), linewidth = 1.1)+
  geom_line(aes(y = q_cournot_wild_own, colour = 'Cournot'), linewidth = 1.1)+
  geom_line(aes(y = q_b_own, colour = 'Bertrand'), linewidth = 1.1)+
  theme_bw()+
  xlab('Stock')+
  ylab("Flow (growth, harvest)")+
  scale_color_manual(values = c('green', 'red', 'grey', 'blue'))

### VII. Rationalizing the way they work for Lemma 2 #####

# Wrong assumption for their proof: 

# They assume that : a_m  = a_w
# These coefficients are computed using the assumption of two producers on the market. 
# Therefore, using a_w would mean : 
# - There is indeed a competitor on the market, and the firm still behaves like a monopoly. 
# In that case, it would mean that the firm behaves like a monopolist on a duopolistic market. Because it 
# operates with residual demand, it produces less. 

# Check : https://www.overleaf.com/5228619579rgmypmgnvbmn

# Therefore, the proof does not compare the right counterfactuals. Let's illustrate that : 

# Define monopoly harvest with direct demand function recomputed
harvest_monop_direct = function(a_m, b_m, c, W, sigma, x){
  y = (a_m - b_m *c)*(sigma^2)*(x^2)/(2*sigma^2*x^2 + 2*b_m*W)
  return(y)
}
a_m = alpha/beta
b_m = 1/beta

# Set data
lemma = data.frame(x, 
                   growth = growth(x,r,k), # Growth
                   monop_aw = harvest_monop_direct(a_w, b_w, c, W, sigma_min, x), # Monopoly using the duopolistic market assumption
                   monop_am = harvest_monop_direct(a_m, b_m, c, W, sigma_min, x), # Monopoly using the monopolistic market assumption
                   q_b_aw = q_b_paper(a_f, a_w, b_f, b_w, c, e, v, W, sigma_min, x)) # Bertrand harvest using the article's formula

colors = c('green', 'blue', 'black')
# Plot
lemma %>% ggplot(aes(x=x))+
  geom_line(aes(y = monop_aw, colour = 'Monopoly on a duopolistic market'), linewidth = 1.1)+
  geom_line(aes(y = q_b_aw, colour = 'Bertrand on a duopolistic market'), linewidth = 1.1)+
  geom_line(aes(y = monop_am, colour = 'Monopoly on a monopolistic market'), linewidth = 1.1, linetype = "dotted")+
  scale_color_manual(values=colors)+
  xlab('x')+
  ylab("Harvest level")+
  theme_bw()
  
# This result shows that if we do not (wrongfully so) assume that a_m = a_w, then their result
# do not hold anymore.
