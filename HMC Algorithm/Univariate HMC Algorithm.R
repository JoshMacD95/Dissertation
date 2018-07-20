#'
#' Further Topics in MCMC
#' 1 Dimension MCMC Algorithm with Hamiltonian Dynamics
#' 
#'

# ==== Preamble

library(coda)
library(mvnfast)
source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')

# ==== Hamiltonian Dynamics Equations


# U is negative log posterior (Posterior is Tareget Distribution)
# grad.U is the derivative of U w.r.t position parameters
# K is function of energy dependent on momentum parameters
# Specific forms of the functions mentioned above should be defined in the 'Functions for HMC.R' file

Univariate.HMC = function(x0, m, L, obs.time, target, K, method=leapfrog, no.its, burn.in){
  # Initial Potenial Energy (-log(posterior))
  time.start = Sys.time()
  U0 = U(x0, target)
  
  x.curr = x0
  U.curr = U0
  
  no.accept = 0
  draws = matrix(nrow = no.its, ncol = 2, byrow = FALSE, data = 0)
  

  for(i in 1:no.its){
    # Momentum Step
    # Draw new momentum value
    rho.curr = rnorm(1, mean = 0, sd = m)
    K.curr = K(rho.curr, m) # Kinetic Energy
    
    # Hamiltonian Dynamics (Leapfrog)
    ham.dyn = numerical.method(x.curr, rho.curr, m, L, obs.time, target, K, method, final = TRUE)
    
    x.prop = ham.dyn[,1]
    rho.prop = -ham.dyn[,2]
    
    U.prop = U(x.prop, target)
    K.prop = K(rho.prop, m)

    # accept-reject step
    if(runif(1) < exp(- U.prop - K.prop + U.curr + K.curr)){
      # Update position and -log posterior
      x.curr = x.prop
      U.curr = U.prop
      rho.curr = rho.prop # Not completely neccesarry as new value drawn
    
      no.accept = no.accept + 1
    }

    # Store New values
    draws[i,1] = x.curr
    draws[i,2] = U.curr
  }
  
  draws = draws[-(1:burn.in),]
  time.finish = Sys.time()
  time.taken = time.finish - time.start
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1])
  scaled.ESS = ESS/L
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(target, K, ESS, accept.rate, x0 = x0, no.its, L, obs.time)
  
  return(list(sample = draws[,1], log.target = -draws[,2], accept.rate = accept.rate, ESS = ESS, scaled.ESS = scaled.ESS))
}

# Simulate Random Normal Data to compare with
#x = rnorm(100000, 0, 1)
#x = sort(x)

# Information about effective sample size
# https://golem.ph.utexas.edu/category/2014/12/effective_sample_size.html


# If obs.time is an even multiple of pi, markov chain will have high, POSITIVE, correlation as 
# proposed values will be close to the current value by hamiltonian dynamics (Full cycles)

# If Obs.time is an odd multiple of pi, markov chain will have high correlation with previous values
# however this correlation alternates between positive and negative 

#  Chosing an obs.time which is between and odd and even multiple of pi will fix these correlation problems
#  giving samples with less correlaation (positive or negative)

#'
#'
#'
#'
#'
#'
#'
# Obs.time too high starts sticking (i.e 100*pi)

# Interesting Result for obs.time = 25*pi
# Explores much faster. Periods of higher variance though
# Seems it occurs for every multiple of 5*pi?

# 15*pi mixes incredibly well, but has moments of sticking.
# Increasing variance of chain

# Problem is Variance parameter, m = 0.5 seems like a good choice. Why?
