### Bayesian analysis with longitudinal studies by indication
### Reagan Mozer and Mark E. Glickman


### This script performs the simulations described in Section 3 of the paper 
### and saves the results.


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")


jags.mod=("jags-model.txt")
jagsscript = 
  cat("
model {
  # priors on parameters
  q1 ~ dunif(0,1);
  q2 ~ dunif(0,1);
  for (j in 1:T.cens){
    p[j] ~ dunif(0,1);
  }
  for (n in 1:N){
  for (t in 1:K){
    prob_ind[n,t] = q1^X[n,t]*q2^(1-X[n,t]);
    Y[n,t] ~ dbern(prob_ind[n,t]);
    time_cross[n,t] = ifelse(Y[n,t]>0, t, T.cens);
  }
  T.obs[n] ~ dsum(min(time_cross[n,]));
  }
  
  for (n in 1:N){
     M_pred[n] = ifelse(T.obs[n]>K, 1, 0); # indicator for ineligibility (i.e., T>K)
     pi[n] = p[T.obs[n]]*(1-M_pred[n]);
     Z.obs[n]~dbern(pi[n]);
  }
  n.control = sum((1-M_pred)*(1-Z.obs));
  n.ineligible = sum(M_pred);
}
",file=jags.mod)


library(rjags)
library(runjags)
library(parallel)
library(plyr)
library(tidyverse)
source("Scripts/00_sim_utils.R")



sim_params=expand.grid(n=c(100,200,500,1000),pX=c(0.1,0.3,0.6,0.9)) %>% arrange(n,pX)

sim_params$seed = c(137550, 30839, 70364, 10217,
                    855, 70568, 58148, 12048,
                    58002, 25680, 44763, 209506,
                    45229, 63265, 139873, 49872)


params = list(q1=0.5, q2=0.25, p=c(0.4,0.5,0.6))


m1 = sim_mod(n=sim_params$n[1],pX=sim_params$pX[1],seed=sim_params$seed[1], 
             mod="jags-model.txt", params=params)

out=purrr::pmap(sim_params, sim_mod, mod="sim1.txt", params=params, save="Results/")
