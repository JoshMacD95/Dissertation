#'
#' Further Topics in MCMC
#' Multivariate HMC Algortihm
#'

# ==== Preamble ====
library(coda)
library(mvnfast)
source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')


# ==== Multivariate Hamiltonian Algorithm ====

Multivariate.HMC = function(x0, m, L, obs.time, 
                            target, K, method = leapfrog, 
                            no.its, burn.in){
  
  time.start = Sys.time()
  
  d = length(x0)
  
  if(length(m) == 1){
    m = rep(m, d)
  }
  if(length(m) != d){
    return(print("Mass vector not same dimension as Position vector"))
  }
  
  # Initial Potential Energy
  U0 = U(x0, target)
  
  x.curr = x0
  U.curr = U0
  
  no.accept = 0
  draws = matrix(data = 0, nrow = no.its + 1, ncol = d + 2)
  
  draws[1,] = c(x.curr, U.curr, U.curr)
  
  
  for(i in 2:(no.its+1)){
    
    # Momentum Step
    rho.curr = rmvn(1, mu = rep(0,d), sigma = diag(m))
    K.curr = K(rho.curr, m)
    H.curr = U.curr + K.curr
    # Hamiltonian Dynamics 
    ham.dyn = numerical.method(x.curr, rho.curr, m, L, obs.time, target, K, method, final = TRUE)
    x.prop = ham.dyn[,1]
    rho.prop = - ham.dyn[,2]
    
    U.prop = U(x.prop, target)
    K.prop = K(rho.prop, m)
    H.prop = U.prop + K.prop
    
    # Accept-recject Step
    if(runif(1) < exp(-H.prop + H.curr)){
      
      # Update position and -log(posteriors)
      x.curr = x.prop
      U.curr = U.prop
      rho.curr = rho.prop
      
      no.accept = no.accept + 1
    }
    
    draws[i,] = c(x.curr, U.curr, H.curr)
  }
  draws = draws[-(1:burn.in),]
  time.finish = Sys.time()
  time.taken = time.finish - time.start
  
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1:d])
  
  # Calculated Scaled ESS which accounts for computation time
  # Divide by number of leapfrog steps L instead of time taken so
  # result is not dependent on computational power.
  scaled.ESS = ESS/L
  
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(target, K, ESS, accept.rate, x0, no.its, L, obs.time)
  
  return(list(sample = draws[,1:d], log.target = -draws[,d+1], hamiltonian = draws[,d+2],  
              ESS = min(ESS), scaled.ESS = min(scaled.ESS), accept.rate = accept.rate))
}
