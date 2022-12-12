### Bayesian analysis with longitudinal studies by indication
### Reagan Mozer and Mark E. Glickman


### This script defines several back-end functions required to run the simulation study described in Section 3 of the paper.

### Problem set up and notation:
##  Study-level (fixed) values
  # N: sample size
  # K: number of study periods
  # pi_X: p(X_t=1) for binary time-varying covariate

## True parameter values
  # q1: P(T=t|X=1)
  # q2: P(T=t|X=0)
  # p: P(Z=1|T=t) for T=1,2,3 



### Data-generating function
generate_data <- function(N, K, pi_X, params){
  
  q1 = params$q1
  q2 = params$q2
  p = params$p

  
  T.cens=K+1  # placeholder value for censored indication times
  
  # 1. Generate a single binary time-varying covariate
  X = matrix(NA, nrow=N, ncol=K)
  X = apply(X,c(1,2),function(x)rbinom(1,1,pi_X))
  
  # 2. Calculate P(T|X) 
  probT_X = q1^X * q2^(1-X)
  
  # 3. Generate T ~ P(T|X)
  Y = apply(probT_X, c(1,2), function(x) rbinom(1,1,x))
  T.ind = get_times(Y)
  
  # 4. Draw Z ~ p(Z|T) with constraint p(Z=1|T>K)=0
  pi_Z = ifelse(T.ind<=K, p[T.ind], 0)
  Z.obs = sapply(pi_Z, function(x) rbinom(1, size=1, prob=x))
  
  Z = ifelse(T.ind==T.cens, NA, Z.obs) # fill in missing assignment for ineligible controls
  
  # 5. Filter observed values based on T and Z
  T.obs=ifelse(Z==1, T.ind, NA)
  Y.obs = filter_Yobs(Y, T.obs)

  
  
  
  # 6. Group summaries
  group = factor(
    dplyr::case_when(Z.obs==1 ~ 1,
              Z.obs==0 & T.ind<=K ~ 2,
              T.ind>K ~ 3),
    labels=c("Treatment","True control", "Ineligible control"))
  
  M = ifelse(T.ind<=K, 0, 1)
  #N.C = sum(group=="True control")
  out = list(dat = list(N=N, K=K, T.cens=K+1, 
                        X=X, Y=Y.obs,
                        T.obs=T.obs, Z.obs=Z.obs),
             T.ind=T.ind, Z=Z, group=group)
  return(out)
}


### Intermediate function to return time of treatment indication
get_times = function(Y){
  N=nrow(Y)
  K=ncol(Y)
  T.cens=K+1
  time.cross=matrix(NA, nrow=N, ncol=K)
  T.ind=c()
  for (n in 1:N){
    for (t in 1:K){
      time.cross[n,t]=ifelse(Y[n,t]==1,t,T.cens)
    }
  }
  T.ind = apply(time.cross,1,min)
  return(T.ind)
}



### Intermediate function to process observed data
filter_Yobs=function(Y, T.ind){
  Y.obs=Y
  Y.obs[is.na(T.ind)]=NA
  N=nrow(Y)
  K=ncol(Y)
  
  for (n in 1:N){
    for (t in 1:K){
      if (!is.na(T.ind[n]) & t>T.ind[n]) Y.obs[n,t]=NA
    }
  }
  return(Y.obs)
}


# Function to randomly impute initial Y values for potential controls and treated when t>T
make_inits=function(T.ind, N, K, imp=T){
  Y=matrix(NA, nrow=N, ncol=K)
  for (n in 1:N){
    for (t in 1:K){
      if (!is.na(T.ind[n]) & t<T.ind[n]){
        Y[n,t]=0
      }
      else if (!is.na(T.ind[n]) & t==T.ind[n]){
        Y[n,t]=1
      }
      
      else if (!is.na(T.ind[n]) & t>T.ind[n] & imp==T){
        Y[n,t]=rbinom(1, 1, 0.5)
      }
      else if (is.na(T.ind[n])){
        Y[n,t]=rbinom(1,1,0.5)
      }
    }
  }
  return(Y)
}

# Randomly draw initial values for each chain
init_mod = function(n=1){
  out=list()
  for (j in 1:n){
    q = runif(2, 0, 1)
    p = runif(4,0,1)
    inits = list(q1=q[1], q2=q[2], p=p)
    out[[j]]=inits
  }
  return(out)
}


# Wrapper function for running multi-factor simulation study with 1,000 iterations
sim_mod = function(n, pX, params, mod,
                   seed=NULL, n.sim=1000, 
                   save=T, dir.out="Results/"){
  
  if(is.null(seed)){seed=Sys.getpid()}
  set.seed(seed)
  
  make_dat = function(){
    d = generate_data(N=n, K=3, pi_X=pX, params=params)
    return(d$dat)
  }
  
  
  m1 = run.jags.study(simulations=n.sim, model=mod, datafunction = make_dat, 
                      n.chains=2,  inits=init_mod(2),
                      targets=params)
  
  
  if (save){
    fname.out = paste0(dir.out, "n",n, "_pX", pX, "_", seed, ".RData")
    save(m1, n, pX, seed, file=fname.out)
  }
  
  return(m1)
  
}
