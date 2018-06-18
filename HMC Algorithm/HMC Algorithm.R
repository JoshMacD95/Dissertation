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

Univariate.HMC = function(q0, V, m, L, obs.time, U, grad.U, K, method, no.its){
  # Initial Potenial Energy (-log(posterior))
  U0 = U(q0, V)
  
  q.curr = q0
  U.curr = U0
  
  no.accept = 0
  draws = matrix(nrow = no.its, ncol = 2, byrow = FALSE, data = 0)
  

  for(i in 1:no.its){
    # Momentum Step
    # Draw new momentum value
    rho.curr = rnorm(1, mean = 0, sd = m)
    K.curr = K(rho.curr, m) # Kinetic Energy
    
    # Hamiltonian Dynamics (Leapfrog)
    ham.dyn = numerical.method(q.curr, rho.curr, m, L, obs.time, grad.U, method, final = TRUE)
    
    q.prop = ham.dyn[,1]
    rho.prop = -ham.dyn[,2]
    
    U.prop = U(q.prop, V)
    K.prop = K(rho.prop, m)
    
    # accept-reject step
    alpha = exp(- U.prop - K.prop + U.curr + K.curr)
    
    if(runif(1) < alpha){
      # Update position and -log posterior
      q.curr = q.prop
      U.curr = U(q.prop, V)
      
      # Store New values
      draws[i,1] = q.curr
      draws[i,2] = U.curr
    
      rho.curr = rho.prop # Not completely neccesarry as new value drawn
      
      no.accept = no.accept + 1
    }
    else{
     draws[i,1] = q.curr
     draws[i,2] = U.curr
    }
  }
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1])
  
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(U, K, ESS, accept.rate, q0 = q0, no.its)
  
  return(draws)
}

test.HMC = Univariate.HMC(q0 = 0, V = 1, m = 1, L = 15.7, obs.time = (3/2)*pi, U = Gaussian, grad.U = grad.Gaussian, K = squared.kinetic, method = leapfrog, no.its = 10000)

L = (3/2)*pi/0.3

# Simulate Random Normal Data to compare with
x = rnorm(100000, 0, 1)
x = sort(x)


# ==== Diagnostics ====

par(mfrow = c(1,2))

plot(test.HMC[,1], type = 'l') # Check for mixing and convergence

acf(test.HMC[,1]) # Check dependence in chain

qqnorm(test.HMC[,1]) # If Data is normal, use this to check distribution

abline(a = 0, b = 1) # The qqplot should follow this line if data is normal

# Plots Histogram of samples and compares with the density of a normal dist.n       
hist(test.HMC[,1], freq = FALSE)
pdf.x = dnorm(x)
lines(x, pdf.x, col = 'blue')       

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


Gaussian(q = c(0,0), V = diag(c(1,1)))

?dmvn

q0 = 0

U0 = Gaussian(q0, V = 1)

draws = matrix(nrow  = 10000, ncol =2, byrow =FALSE, data = 0)
