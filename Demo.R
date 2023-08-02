library(psych)
library(pscl)
library(dplyr)
library(mvtnorm)
library(foreach)
library(doParallel)
library(ranger)
library(tidyverse)
library(BNSP)

# The data shall be scaled as 0 -- J-1 if it is a J-point scale
y = read.csv('Data_BFI.csv') %>% filter(Condition==3) %>% select(!all_of('Condition'))
# number of items in each dimension
## The example is five dimensions containing 10 items, respectively
pois = rep(10,5)
# the corresponding dimension for each item
## The example is items 1-9 are in dimension 1; items 10-19 are in dimension 2; items 21-29 are in dimension 3; items 31-39 are in dimension 4; items 41-49 are in dimension 5.
t_pois = rep(1:5,each=10)

# J: J-point scale
# nchain: number of chains, must be larger than 1
# iter: total number of iterations
# burnin: number of burn-in
# thin: number of thining
T1 = Sys.time()
mod = RDT(y,J=5,pois=pois,t_pois=t_pois,nchain=3,iters=21000,burnin=7000,thin=14)
Sys.time() - T1

# ReTrain: update the earlier model if its parameters did not converge; default is F
# If ReTrain is T, ``mod'' must be the result of the earlier model
T2 = Sys.time()
mod = RDT(y,J=5,pois=pois,t_pois=t_pois,nchain=3,iters=15000,burnin=5000,thin=10,ReTrain=T,mod=mod)
Sys.time() - T2

# EAP estimates
## Alpha, Gamma, Delta, Omega, Tau, ThetaR, ThetaD, ThetaT, Sigma11, Sigma22
mod$Estimate$Alpha

# R-Square
## Alpha, Gamma, Delta, Omega, Tau, ThetaR, ThetaD, ThetaT
mod$RSquare$Delta

# MCMC plot
## Alpha, Gamma, Delta, Omega, Tau, ThetaD, ThetaT
### i: item or individual
plotRDTMCMC(mod, x = 'Omega', i=1)
## ThetaR
### d: dimension
plotRDTMCMC_M(mod, x = 'ThetaR', d=1, i=1)

# DIC value: LogL: -2logLikelihood; DIC: DIC
DIC(y,mod,nchain=3,t_pois=t_pois)
