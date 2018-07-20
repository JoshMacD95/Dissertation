#'
#' Further Topics in MCMC
#' Random Walk Metropolis Algorithm (Symmetric Proposal)
#'
#'

# ==== Preamble ====
library(coda)
source('Other MCMC methods/Functions for MCMC.R')

# ==== RWM Algorithm ====
# A Random Walk Metropolis Algorithm for Symmetric proposals
RWM = function(x0, target, proposal, no.its, lambda){
  
  # Record starting time
  time.start = system.time()
  
  # Store Dimension of problem
  d = length(x0)
  
  # Assign the intial value as the current value in the Markov chain
  x.curr = x0
  
  # Calculate the current values density value
  target.curr = target(x0, dist = target)
  
  # Track the number of propsals accepted
  accept = 0
  
  # Create matrix to store samples and density values
  draws = matrix(0, nrow = no.its + 1, ncol = d + 1)
  samples[1,] = c(x.curr, target.curr)
  
  for(i in 2:no.it+1){
    # Propose new value using chosen symmetric proposal 
    x.prop = x.curr + lambda*proposal(dist = proposal)
    target.prop = density.fn(x = x.prop, dist = target)
    
    # Accept/Reject Step with Metropolis-Hastings Acceptance probability
    if(runif(1) < target.prop/target.curr){
      # Set proposed value as current value in Markov chain if accepted
      x.curr = x.prop
      target.curr = target.prop
      accept = accept + 1
    }
    # Otherwise, the new current value is the same as the previous current value.
    
    # Store current value and density value
    samples[i,] = c(x.curr, target.curr)
  }
  time.taken = system.time - time.start
  ESS = effective.size(samples[,1:d])
  accept.rate = accept/no.its
  
  output.RWM(target, proposal, lambda, ESS, accept.rate, x0)
  

}