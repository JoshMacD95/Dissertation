#'
#' Further Topics in MCMC
#' Posterior Functions for HMC
#'

# ==== Potential Energy Function ==== 
# (-log(target))

Gaussian = function(q, sigma = 1){
  return(-log(dnorm(q,0,sigma)))
}

grad.Gaussian= function(q){
  return(q)
}

# Kinetic Energy Function
# Acts as a proposal in HMC
squared.kinetic = function(rho,m){
  return(0.5*rho^2/m)
} 

# ==== Output for HMC Algorithm ====

output = function(U, K, ESS, accept.rate, q0, no.its){
  # Problem with Target and Kinetic function names
  print(paste(c("Target: ", U), sep = "", collapse = ""))
  print(paste(c("Kinetic Energy Function:", K)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size:", ESS)), sep = "", collapse = "")
  print(paste(c("Acceptance Rate:", accept.rate)), sep = "", collapse = "")
  print(paste(c("Initial Conditions:", q0)), sep = "", collapse = "")
  print(paste(c("Number of Iterations:", no.its)), sep = "", collapse = "")
}


# Total Energy (This should stay constant by the conservation of energy)
#H = function(q, rho, m){
#  return(U(q) + K(rho,m))
#}
