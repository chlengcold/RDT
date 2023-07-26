library(psych)
library(pscl)
library(dplyr)
library(mvtnorm)
library(foreach)
library(doParallel)
library(ranger)
library(tidyverse)

# The data shall be scaled as 0 -- J-1 in a J-point item.
y = read.csv('Data_BFI.csv') %>% filter(Condition==3) %>% select(!all_of('Condition'))
# number of scale
J= 5 
# number of items in each dimension
pois = rep(10,5)
# the corresponding dimension for each item
t_pois = rep(1:5,each=10)

mod = RDT(y,J,pois=pois,t_pois=t_pois,nchain=3,iters=21000,burnin=7000,thin=14)

# EAP estimates
## Alpha, Gamma, Delta, Omega, Tau, ThetaR, ThetaD, ThetaT, Sigma11, Sigma22
mod$Estimate$Alpha

# R-Square
## Alpha, Gamma, Delta, Omega, Tau, ThetaR, ThetaD, ThetaT
mod$RSquare$Alpha

# MCMC plot
## Alpha, Gamma, Delta, Omega, Tau, ThetaD, ThetaT
### i: item or individual
plotRDTMCMC(mod, x = 'Alpha', i=1)
## ThetaR
### d: dimension
plotRDTMCMC_M(mod, x = 'ThetaR', d=1, i=1)

# DIC value: LogL: -2logLikelihood; DIC: DIC
DIC(y,mod,nchain=3,t_pois=t_pois)
