
#Author: Erendira Aceves
#Goal: Dispersal and Sustainability in Spatial management
#Date: Feb 7

################################# PARAMETERS #########################

#Weight
alpha <- 3 #Fixed parameter of weight function
beta <- .5

#Fecundity
mu <- 3 #Fixed parameter of fecundity function. RAte (slope)
b <- 3 #Fixed parameter of fecunfity function. Intecept

#Price
gamma <-0
d <- .05

# Recruiment
h <- .7#Most common value in marine species (Hilborn)
R0 <-10000 #Initial number of recruits


#Popoulation growth
m <- 0.05 #Natural mortality
amax = 20 #Because this way I am setting the mortality to 1/lifespan
increment = .5
a <- seq(0,amax,by=increment) #Ages
years_max <- 100
R0 <- 1000 #Initial number of Recruits
winf <- 500

## Choice variable: amu

################################# FUNCTIONS #########################

Weight <- function(alpha,beta,a,winf){
  weights <- rep(NA, length(a))
  weights <- winf*((1-exp(-beta*a))^alpha)
  return(weights)
}


Fecundity <- function(mu,b,alpha,beta,a,m,winf){
  w <- Weight(alpha=alpha,beta=beta,a=a,winf=winf)
  amat_loc <- which.max(diff(w))
  fec <- (mu*w)+b
  fec[1:amat_loc-1] <- 0
  SPR_0 <- sum(fec*exp(-m))
  return(list(fec,SPR_0))
}

Price <- function(alpha,beta,a,gamma,winf){
  w <- Weight(alpha,beta,a,winf)
  price <- exp(gamma*w)
  price<- rep(NA,length(a))
  price <- (exp(gamma*w))
  return(price)
}


Recruits <- function(N,R0,h,mu,b,alpha,beta,a,m,increment,winf,...){
  # The fecundity of each age is devided into smaller intervals when the increment of age is <1
  a_rounded <- seq(0,amax,by=1) 
  Fec_temp <- Fecundity(mu=mu,b=b,alpha=alpha,beta=beta,a=a_rounded,m=m,winf=winf)[[1]]
  SPR0 <- Fecundity(mu=mu,b=b,alpha=alpha,beta=beta,a=a_rounded,m=m,winf=winf)[[2]]
  interval <- (1/increment)
  Fec_values <- Fec_temp/interval
  Fec <- rep(0,length(a))
  #### This loop devided the fecundity values to match the number of steps (interval)
  for (i in 1:(length(a_rounded)-1)){
    Fec_loc <- which(a==i)
    Fec[(Fec_loc-interval):(Fec_loc-1)] <- rep(Fec_values[i],interval) 
  }
  betarec <- ((5*h)-1)/(4*h*R0)
  alpharec <- SPR0*((1-h)/(4*h))
  # N (number of indivudal in each age class) and Fec are vectors of the same length
  # The amount of larvae produced by each age class after deensity dependence is:
  Larvae.age <- Fec*N
  Larvae.year <- sum(Larvae.age)
  #The amount of recruits produced every year after density dependence 
  Recruits.Year <- (Larvae.year)/(alpharec+(betarec*Larvae.year))
  return(Recruits.Year)
}

########################### Survaival

Surv <- function(m,amu,a,increment) {
  interval <- 1/increment
  amu_loc <- which(a==amu)
  surv <- rep(NA,length(a))
  mortality <- 1- ((1-m)^(1/interval))
  surv[1:amu_loc-1] <- exp((-mortality/interval))
  surv[amu_loc:length(a)] <- 0
  return(surv)
}


####################################### POPULATION GROWTH

### Intial population

IP <- function(a,R0,amax,m,increment,...){
  amu_Ipop <- amax
  I.pop <- rep(NA,length(a))
  I.pop[1] <- R0
  for (i in 2:length(a)){
    I.pop[i] <- Surv(m=m,amu=amu_Ipop,a=a,increment=increment)[i]*I.pop[i-1]
  }
  return(I.pop)
}


### Population Matrix
Pop <- function (a,R0,amax,m,amu,increment,mu,h,winf,...) {
  M<- matrix(NA,length(a),years_max)
  M[,1] <- IP(a=a,R0=R0,amax=amax,m=m,increment=increment)
  R1 <- Recruits(N=M[,1],R0=R0,h=h,mu=mu,b=b,alpha=alpha,beta=beta,a=a,m=m,increment=increment,winf=winf)
  R <- R1
  for (i in 2:years_max){
    M[1,i] <- R
    for (j in 2:length(a)){
      M[j,i] <- M[j-1,i-1]*Surv(m=m,amu=amu,a=a,increment=increment)[j]
    }
    R<- Recruits(N=M[,i],R0=R0,h=h,mu=mu,b=b,alpha=alpha,beta=beta,a=a,m=m,increment=increment,winf=winf)
  }
  return(M)
}

Harvest <- function(a,R0,amax,m,amu,alpha,beta,increment,mu,h,winf,...){
  amu_loc <- which(a==amu)
  w <- Weight(alpha=alpha,beta=beta,a=a,winf=winf)
  yield_temp <- Pop(a=a,R0=R0,amax=amax,m=m,amu=amu,increment=increment,mu=mu,h=h,winf=winf)[(amu_loc-1),]
  yield <-yield_temp*w[amu_loc-1]
  return(yield)}

############################# NPV
NPV <- function(years_max,d,alpha,beta,a,R0,amax,m,amu,gamma,increment,mu,h,winf,...){
  amu_loc <- which(a==amu)
  p <- Price(alpha=alpha,beta=beta,a=a,gamma=gamma,winf=winf)
  yield <- Harvest(a=a,R0=R0,amax=amax,m=m,amu=amu,alpha=alpha,beta=beta,increment=increment,mu=mu,h=h,winf=winf)
  profits_year <- yield*p[amu_loc-1]
  years <- seq(1,years_max)
  Dis.profits <- profits_year*(1/((1+d)^years))
  NPV <- sum(Dis.profits)
  return(NPV)
}

############################### FINDING THE BEST AGE OF FIRST HARVEST


#Optimal age of first harvest for Social Planner  
amuSP <- function(years_max,d,alpha,beta,R0,amax,m,gamma,increment,mu,h,winf,...){
  a <- seq(0,amax,by=increment)
  NPV_amu <- rep(0,length(a)) 
  for (i in 2:length(a)){
    NPV_amu[i] <- NPV(years_max=years_max,d=d,alpha=alpha,beta=beta,a=a,R0=R0,amax=amax,m=m,amu=a[i],gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)
  }
  maxNPV <- max(NPV_amu)
  optage_loc <- which(NPV_amu==maxNPV)
  opt_ageSP <- a[optage_loc]
  return(opt_ageSP)
}

# Optimal age of first harvest for TURF owner ### Faustman  
amuTO<- function(alpha,beta,d,m,gamma,increment,amax,winf){
  a <- seq(0,amax,by=increment)
  interval <- 1/increment
  profits <- rep(NA,length(a))
  w <- Weight(alpha,beta,a,winf)
  price <- Price(alpha=alpha,beta=beta,a=a,gamma=gamma,winf=winf)
  for (i in 1:length(a)){
    profits[i] <- price[i]*w[i]*exp(-(d+m)*a[i])
  }
  max_profits <- max(profits)
  optage_loc <- which(profits==max_profits)
  opt_ageTO  <- a[optage_loc]
  return(opt_ageTO)
}




################################ 

#Age of first harvest SP
optimal.age.SP <- amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)

#Age of firts harvest TO
optimal.age.TO <- amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)

age.maindiff <- optimal.age.SP-optimal.age.TO 
age.maindiff <- age.maindiff/optimal.age.SP
percentdiff <- age.maindiff*100

### Main diff NPV
NPV.SP <- NPV(years_max=years_max,d=d,alpha=alpha,bet=beta,a=a,R0=R0,amax=amax,m=m,amu=optimal.age.SP,gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)
NPV.TO <- NPV(years_max=years_max,d=d,alpha=alpha,bet=beta,a=a,R0=R0,amax=amax,m=m,amu=optimal.age.TO,gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)

NPV.maindiff <- NPV.SP-NPV.TO
NPV.maindiff <- NPV.maindiff/NPV.SP
percentdiffNPV <- NPV.maindiff*100




######################################## THE EFFECT OF A PRICE PREMIUM

#install.packages("foreach") ### Allows parallel execution
library(foreach)

gamma_values <- seq(0,1,by=.1)
gamma_amuSP <-  rep(NA,length(gamma_values))
gamma_amuTO <- rep(NA,length(gamma_values))

gamma_amuSP <- foreach (i=1:length(gamma_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma_values[i],increment=increment,mu=mu,h=h,winf=winf)[1]
}

gamma_amuTO <- foreach (i=1:length(gamma_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma_values[i],increment=increment,amax=amax,winf=winf)[1]
}

gammadiff <- gamma_amuSP-gamma_amuTO
gammadiff <- (gammadiff/gamma_amuSP)*100

#### NPV with different growth rate
NPV_gamma_amuSP <- rep(NA,length(gamma_values))
NPV_gamma_amuTO <- rep(NA,length(gamma_values))

NPV_gamma_amuSP <- foreach (i=1:length(gamma_values),.combine='c') %do%{
  NPV(years_max=years_max,d=d,alph=alpha,beta=beta,a=a,R0=R0,amax=amax,m=m,amu=gamma_amuSP[i],gamma=gamma_values[i],increment=increment,mu=mu,h=h,winf=winf)
}

NPV_gamma_amuTO <- foreach (i=1:length(gamma_values),.combine='c') %do% {
  NPV(years_max=years_max,d=d,alpha=alpha,beta=beta,a=a,R0=R0,amax=amax,m=m,amu=gamma_amuTO[i],gamma=gamma_values[i],increment=increment,mu=mu,h=h,winf=winf)
}

gammadiffNPV <- NPV_gamma_amuSP-NPV_gamma_amuTO
gammadiffNPV <- (gammadiffNPV/NPV_gamma_amuSP)*100

attach(mtcars)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),widths = c(.8,.8),
       heights = c(.8,.8),respect=F)
par(oma=c(2,4,2,4)) #bottom, left, top, right

plot(gamma_values,gamma_amuSP,type='l',col='turquoise4',ylim=c(0,amax),ylab='Age of First Harvest',xlab='Price Increase Rate',lwd=3,xaxt='n')
lines(gamma_values,gamma_amuTO,col='turquoise',lwd=3,lty=3)
mtext("Age of 1st Harvest",side=2,line=2.5)
legend(.85,15,'1A)',box.col= "white",cex=.8)
legend(.40,10,c('Social Planner','TURF Owner'),box.col="white",lwd=3,lty=c(1,3),col=c('turquoise4','turquoise'),y.intersp=1.5,cex=.8)

plot(gamma_values,gammadiffNPV, type='l',xlab='', col='chocolate4', ylim=c(0,100),ylab='',lwd=3,yaxt='n')
axis(4)
lines(gamma_values, gammadiff, col='chocolate3', main='Age of First Harvest',ylim=c(0,100),lty=3,lwd=3)
mtext("Perecent Difference",side=4,line=2.5)
legend(.85,.99,'1B)',box.col= "white",cex=.8)
legend(.40,.9,c('NPV','Age of 1st Harvest'),box.col="white",lwd=3,lty=c(1,3),col=c('chocolate4','chocolate3'),y.intersp=1.5,cex=.8)
mtext("Price Increase Rate",side=1,line=2.5)






############################### Sensitivity anlysis

################ Varying the individual growth rate
beta_values <- seq(.1,1,by=.1)
beta_amuSP <-  rep(NA,length(beta_values))
beta_amuTO <- rep(NA,length(beta_values))

beta_amuSP <- foreach (i=1:length(beta_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta_values[i],R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,h=h,mu=mu,winf=winf)[1]
}

beta_amuTO <- foreach (i=1:length(beta_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta_values[i],d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

betadiff <- beta_amuSP-beta_amuTO
betadiff <- (betadiff/beta_amuSP)*100


########################## Mortality Rate

m_values <- seq(0.01,.1,by=.01)
m_amuSP <-  rep(NA,length(m_values))
m_amuTO <- rep(NA,length(m_values))

m_amuSP <- foreach (i=1:length(m_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=amax,m=m_values[i],gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)[1]
}

m_amuTO <- foreach (i=1:length(m_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m_values[i],gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

mdiff <- m_amuSP-m_amuTO
mdiff <- (mdiff/m_amuSP)*100


################ Varying the Discount Rate

delta_values <- seq(.01,.1,by=.01)
delta_amuSP <-  rep(NA,length(delta_values))
delta_amuTO <- rep(NA,length(delta_values))

delta_amuSP <- foreach (i=1:length(delta_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=delta_values[i],alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)[1]
}

delta_amuTO <- foreach (i=1:length(delta_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=delta_values[i],m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

deltadiff <- delta_amuSP-delta_amuTO
deltadiff <- (deltadiff/delta_amuSP)*100

############################# fecundity

mu_values <- seq(1,10,by=1)
mu_amuSP <-  rep(NA,length(mu_values))
mu_amuTO <- rep(NA,length(mu_values))

mu_amuSP <- foreach (i=1:length(mu_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,mu=mu_values[i],h=h,winf=winf)[1]
}

mu_amuTO <- foreach (i=1:length(mu_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

mudiff <- mu_amuSP-mu_amuTO
mudiff <- (mudiff/mu_amuSP)*100

################################ Steepness

h_values <- seq(.1,.9,by=.1)
h_amuSP <-  rep(NA,length(h_values))
h_amuTO <- rep(NA,length(h_values))

h_amuSP <- foreach (i=1:length(h_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,mu=mu,h=h_values[i],winf=winf)[1]
}

h_amuTO <- foreach (i=1:length(h_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

hdiff <- h_amuSP-h_amuTO
hdiff <- (hdiff/h_amuSP)*100

################################ amax

amax_values <- seq(1,50,by=5)
amax_amuSP <-  rep(NA,length(amax_values))
amax_amuTO <- rep(NA,length(amax_values))

amax_amuSP <- foreach (i=1:length(amax_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,amax=50,m=m,gamma=gamma,increment=increment,mu=mu,h=h,winf=winf)[1]
}

amax_amuTO <- foreach (i=1:length(amax_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

amaxdiff <- amax_amuSP-amax_amuTO
amaxdiff <- (amaxdiff/amax_amuSP)*100

############################### Winf

winf_values <- seq(100,1000,by=100)
winf_amuSP <-  rep(NA,length(winf_values))
winf_amuTO <- rep(NA,length(winf_values))

winf_amuSP <- foreach (i=1:length(winf_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0,m=m,gamma=gamma,increment=increment,mu=mu,h=h,amax=amax,winf=winf_values[i])[1]
}

winf_amuTO <- foreach (i=1:length(winf_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf_values[i])[1]
}

winfdiff <- winf_amuSP-winf_amuTO
winfdiff <- (winfdiff/winf_amuSP)*100


##################################

alpha_values <- seq(1,10,by=1)
alpha_amuSP <-  rep(NA,length(alpha_values))
alpha_amuTO <- rep(NA,length(alpha_values))

alpha_amuSP <- foreach (i=1:length(alpha_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha_values[i],beta=beta,R0=R0,m=m,gamma=gamma,increment=increment,mu=mu,h=h,amax=amax,winf=winf)[1]
}

alpha_amuTO <- foreach (i=1:length(alpha_values),.combine='c') %do% {
  amuTO(alpha=alpha_values[i],beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

alphadiff <- alpha_amuSP-alpha_amuTO
alphadiff <- (alphadiff/alpha_amuSP)*100

###########################
R0_values <- seq(100,5000,by=500)
R0_amuSP <-  rep(NA,length(R0_values))
R0_amuTO <- rep(NA,length(R0_values))

R0_amuSP <- foreach (i=1:length(R0_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=d,alpha=alpha,beta=beta,R0=R0_values[i],m=m,gamma=gamma,increment=increment,mu=mu,h=h,amax=amax,winf=winf)[1]
}

R0_amuTO <- foreach (i=1:length(R0_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=d,m=m,gamma=gamma,increment=increment,amax=amax,winf=winf)[1]
}

R0diff <- R0_amuSP-R0_amuTO
R0diff <- (R0diff/R0_amuSP)*100



############################################## SENSITIVITY Plots

#Sensitivity to rates
#attach(mtcars)
layout(matrix(c(1,2,3,4,5,6,7,8,9,10), 3, 3, byrow = TRUE),widths = c(.8,.8,.8),
       heights = c(1,1,1),respect=F)
par(oma=c(2,2,2,2)) #bottom, left, top, right

par(mar=c(5,4,0,0))
plot(beta_values, betadiff, type='l',xlab='Individual Growth Rate',ylab='Percent Difference',lwd=2,lty=3,ylim=c(0,100))
legend('topright','A)',box.col= "white")

par(mar=c(5,4,0,0))
plot(m_values, mdiff, type='l',xlab='Mortality',ylim=c(0,100),lwd=2,yaxt='n',xlim=c(.01,.1),lty=3,ylab='')
legend('topright','B)',box.col= "white")

par(mar=c(5,4,0,0))
plot(delta_values,deltadiff, type='l',xlab='Discount Rate', ylim=c(0,100),lwd=2,yaxt='n',lty=3,ylab='')
legend('topright','C)',box.col= "white")

par(mar=c(4,4,1,0))
plot(mu_values, mudiff, type='l',xlab='Fecundity Rate', ylim=c(0,100),ylab='Percent Difference',lwd=2,lty=3)
legend('topright','D)',box.col= "white")

par(mar=c(4,4,1,0))
plot(h_values,hdiff, type='l',xlab='Steepness',ylim=c(0,100),lwd=2,yaxt='n',lty=3,ylab='')
legend('topright','E)',box.col= "white")

par(mar=c(4,4,1,0))
plot(amax_values,amaxdiff, type='l',xlab='Longevity', ylim=c(0,100),lwd=2,yaxt='n',lty=3,ylab='')
legend('topright','F)',box.col= "white")


par(mar=c(4,4,1,0))
plot(winf_values, winfdiff, type='l',xlab='Maximum weight', ylim=c(0,100),ylab='Percent Difference',lwd=2,lty=3)
legend('topright','G)',box.col= "white")

par(mar=c(4,4,1,0))
plot(alpha_values,alphadiff, type='l',xlab='Intercep of Fecundity',ylim=c(0,100),lwd=2,yaxt='n',lty=3,ylab='')
legend('topright','H)',box.col= "white")

par(mar=c(4,4,1,0))
plot(R0_values,R0diff, type='l',xlab='Intial Recruits', ylim=c(0,100),lwd=2,yaxt='n',lty=3,ylab='')
legend('topright','I)',box.col= "white")




################################# LOCO
m <- 0.1 #Natural mortality
amax = 10 #maximum age
amat=4.5 #En realidad es 4.5 pero lo redondie porque la funcion de reclutas no reconoce el 0.5 (usa a_rounded) 
kvalues <- c(.161,.203,.215,.347,.362,.12,.38,.52,.17,.209,.32,.55)
loco_values <- seq(0,.1,0.01) #Discount rate values
#Ref IMARPE 1996, RAVI MARAVI 1997
beta<- median(kvalues)
winf <- 250
a<- seq(0,amax,.5)

Price <- function(alpha,beta,a,gamma,winf){
  w <- Weight(alpha=alpha,beta=beta,a=a,winf=winf)
  price <- -8.662e-01+(1.155e-01*w)+(-4.389e-04*w^2)+(5.386e-07*w^3)
  return(price)
}

############### For all the calculated growth rates ofloc0
loco_amuSP <-  rep(NA,length(loco_values))
loco_amuTO <- rep(NA,length(loco_values))

loco_amuSP <- foreach (i=1:length(loco_values),.combine='c') %do% {
  amuSP(years_max=years_max,d=loco_values[i],alpha=alpha,beta=beta,R0=R0,amax=amax,m=m,gamma=gamma,increment=increment,amat=amat,h=h,mu=mu,winf=winf)[1]
}

loco_amuTO <- foreach (i=1:length(loco_values),.combine='c') %do% {
  amuTO(alpha=alpha,beta=beta,d=loco_values[i],m=m,gamma=gamma,increment=increment,winf=winf,amax=amax)[1]
}

locodiff <- loco_amuSP-loco_amuTO
locodiff <- locodiff/loco_amuSP

#### NPV with different growth rate
NPV_loco_amuSP <- rep(NA,length(loco_values))
NPV_loco_amuTO <- rep(NA,length(loco_values))

NPV_loco_amuSP<- foreach (i=1:length(loco_values),.combine='c') %do% {
  NPV(years_max=years_max,d=loco_values[i],alpha=alpha,beta=beta,a=a,R0=R0,amax=amax,m=m,amu=loco_amuSP[i],gamma=gamma,increment=increment,amat=amat,h=h,mu=mu,winf=winf)
}

NPV_loco_amuTO<- foreach (i=1:length(loco_values),.combine='c') %do% {
  NPV(years_max=years_max,d=loco_values[i],alpha=alpha,beta=beta,a=a,R0=R0,amax=amax,m=m,amu=loco_amuTO[i],gamma=gamma,increment=increment,amat=amat,h=h,mu=mu,winf=winf)
}


locodiffNPV <- NPV_loco_amuSP-NPV_loco_amuTO
locodiffNPV <- locodiffNPV/NPV_loco_amuSP

dev.off()
plot(loco_values,locodiffNPV, type='l', col='darkolivegreen4', ylim=c(0,1),ylab='Perecent Difference',xlab='Discount Rate',lwd=3)
lines(loco_values, locodiff, col='darkolivegreen3', main='Age of First Harvest',ylim=c(0,1),lty=3,lwd=3)
legend(.06,.99,c('NPV','Age of 1st Harvest'),col=c('darkolivegreen4','darkolivegreen3'),lty=c(1,3),lwd=3,box.col= "white")



