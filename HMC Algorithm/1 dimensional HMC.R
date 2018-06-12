#'
#' Further Topics in MCMC
#' MCMC Algorithm with Hamiltonian Dynamics
#' 1-dimensional HMC 
#'

# ==== Preamble

library(coda)
setwd('H:/MSc Statistics/Dissertation/Hamiltonian Monte Carlo/1 Dimensional HMC Algorithm')

source('Functions for HMC.R')
# ==== Hamiltonian Dynamics Equations


# U is negative log posterior (Posterior is Tareget Distribution)
# grad.U is the derivative of U w.r.t position parameters
# K is function of energy dependent on momentum parameters
# Specific forms of the functions mentioned above should be defined in the 'Functions for HMC.R' file

HMC = function(U, grad.U, K, timestep, obs.time, q0, no.its, m){
  
  q.curr = q0
  no.accept = 0
  draws = c()
  for(i in 1:no.its){
    
    # Momentum Step
    # Draw new momentum value
    
    rho.curr = rnorm(1, 0, m)
    
    # Hamiltonian Dynamics (Leapfrog)
    q = q.curr
    rho = rho.curr
    
    U.curr = U(q.curr)
    K.curr = K(rho.curr, m)
    
    # Leapfrog 
    j = 1
    while(j*timestep < obs.time){
      rho_halfeps = rho - 0.5*timestep*grad.U(q)
      q = q + timestep*rho_halfeps/m
      rho = rho_halfeps - 0.5*timestep*grad.U(q)
      j = j + 1
    }
    
    q.prop = q
    rho.prop = -rho
    
    U.prop = U(q.prop)
    K.prop = K(rho.prop, m)
    
    
    # accept-reject step
    alpha = exp(- U.prop - K.prop + U.curr + K.curr)
    
    if(runif(1) < alpha){
      q.curr = q.prop
      draws[i] = q.curr
      # Not completely neccesarry as new value drawn
      rho.curr = rho.prop
      no.accept = no.accept + 1
    }
    else{
     draws[i] = q.curr
    }
  }
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws)
  
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(U, K, ESS, accept.rate, q0 = q0, no.its)
  
  return(draws)
}

test.HMC = HMC(U = Gaussian, grad.U = grad.Gaussian, K = squared.kinetic, timestep = 0.1, obs.time = 2*pi, q0 = 0, no.its = 10000, m = 0.6)


x = rnorm(100000, 0, 1)
x = sort(x)
qqnorm(x)
par(mfrow = c(1,1))
abline(a = 0, b = 1)

plot(test.HMC, type = 'l')
acf(test.HMC)
qqnorm(test.HMC)
abline(a = 0, b = 1)
       
hist(test.HMC, freq = FALSE)
pdf.x = dnorm(x)
lines(x, pdf.x, col = 'blue')       
# https://golem.ph.utexas.edu/category/2014/12/effective_sample_size.html
)

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





