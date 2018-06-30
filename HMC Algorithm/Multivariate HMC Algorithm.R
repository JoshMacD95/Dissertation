#'
#' Further Topics in MCMC
#' Multivariate HMC Algortihm
#'

# ==== Preamble ====
library(coda)
library(mvnfast)
source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')

# ==== Multivariate Hamiltonian Algorithm ====

Multivariate.HMC = function(q0, m, L, obs.time, 
                            U, grad.U, K, method = leapfrog, 
                            no.its, burn.in, output.type = "all"){
  
  time.start = Sys.time()
  
  d = length(q0)
  
  if(length(m) == 1){
    m = rep(m, d)
  }
  if(length(m) != d){
    return(print("Mass vector not same dimension as Position vector"))
  }
  
  # Initial Potential Energy
  U0 = U(q0)
  
  q.curr = q0
  U.curr = U0
  
  no.accept = 0
  draws = matrix(data = 0, nrow = no.its + 1, ncol = d + 2)
  
  draws[1,] = c(q.curr, U.curr, U.curr)
  
  
  for(i in 2:(no.its+1)){
    
    # Momentum Step
    rho.curr = rmvn(1, mu = rep(0,d), sigma = diag(m))
    K.curr = K(rho.curr, m)
    H.curr = U.curr + K.curr
    # Hamiltonian Dynamics 
    ham.dyn = numerical.method(q.curr, rho.curr, m, L, obs.time, grad.U, method, final = TRUE)
    q.prop = ham.dyn[,1]
    rho.prop = - ham.dyn[,2]
    
    U.prop = U(q.prop)
    K.prop = K(rho.prop, m)
    H.prop = U.prop + K.prop
    
    # Accept-recject Step
    if(runif(1) < exp(-H.prop + H.curr)){
      
      # Update position and -log(posteriors)
      q.curr = q.prop
      U.curr = U.prop
      rho.curr = rho.prop
      
      no.accept = no.accept + 1
    }
    
    draws[i,] = c(q.curr, U.curr, H.curr)
  }
  draws = draws[-(1:burn.in),]
  time.finish = Sys.time()
  time.taken = time.finish - time.start
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1:d])
  ESS.per.sec = ESS/as.numeric(time.taken)
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(U, K, ESS, accept.rate, q0, no.its)
  if(output.type == "all"){
    return(list(sample = draws[,1:d], log.target = -draws[,d+1], hamiltonian = draws[,d+2],  
                ESS = ESS, ESS.per.sec = ESS.per.sec))
  }
  if(output.type == "ESS"){
    return(ESS)
  }
  if(output.type == "ESS/sec"){
    return(ESS.per.sec)
  }
}
