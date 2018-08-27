#'
#'
#' RWM Algorithm 
#'
#'

# ==== Preamble ====

source("RWM Algorithm/Functions for RWM.R")
library(coda)

# ==== The Algorithm ====
Multivariate.RWM=function(no.its, lambda,target="Std.Gaussian",proposal="Std.Gaussian", x0, prop.V, burn.in){
  
  time.start = Sys.time()
  d = length(x0)
  draws = matrix(data = 0, nrow = no.its + 1, ncol = d + 1)
  no.accept=0
  
  x.curr=x0
  log.density.curr=log.density(x.curr,target)

  draws[1,] = c(x.curr, log.density.curr)
  for(i in 1:no.its){
    z=propose(proposal, d, V = prop.V)
    # Proposed Value is Current Value plus random noise from proposal dist.n
    # Lambda controls the 
    x.prop = x.curr + lambda*z
    
    log.density.prop = log.density(x.prop,target)
    log.alpha = log.density.prop - log.density.curr
    
    u = runif(1)
    
    if(log(u) <= log.alpha){ # accept
      x.curr = x.prop
      log.density.curr = log.density.prop
      no.accept = no.accept+1
    }
    
    draws[i,1:d] = x.curr # always store
  }
  
  draws = draws[-(1:burn.in),]
  
  time.finish = Sys.time()
  time.taken = time.finish - time.start
  
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1:d])
  
  # Calculated Scaled ESS which accounts for computation time
  # Divide by number of leapfrog steps L instead of time taken so
  # result is not dependent on computational power.
  scaled.ESS = ESS/as.numeric(time.finish - time.start)
  
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  print("1D RWM")
  print("======")
  
  output.RWM(target, proposal, lambda, ESS, accept.rate, no.its)
  return(list(sample = draws[,1:d], target.density = exp(draws[,d+1]),  
              ESS = min(ESS), scaled.ESS = min(scaled.ESS), accept.rate = accept.rate))
}






       