#'
#' Further Topics in MCMC
#' Functions for Univariate HMC
#'

# ==== Preamble ====
library(mvnfast)

# ==== Potential Energy Function ==== 
# (-log(target))

# == Gaussian ==
Gaussian = function(q){
  return(-dnorm(q, 0, sd = 1, log = TRUE))
}

grad.Gaussian = function(q){
  return(q)
}

# == Logistic ==

Logistic = function(q){
  return(-log(exp(q)/(1+exp(q))^2))
}

grad.Logistic = function(q){
  return(-1 + 2*exp(q)/(1+exp(q)))
}

# == Multivariate Gaussian ==
"
U = function(x, mu, V, target){
  if(target = "Gaussian"){
    return(-dmvn(x, mu, sigma = V, 
                 log = TRUE))
  }
  else if(target = "Logistic"){
    return(sum(-log(exp(x)/(1+exp(x))^2)))
  }
  else{
    stop("Please give a valid target:
         'Gaussian'
         'Logistic'
         ")
  }
}
"
"
grad.U = function(x, mu, V, target){
  
  if(target = "Gaussian"){
    return(x%*%solve(V))
  }
  else if(target = "Logistic"){
    return(-1 + 2*exp(x)/(1+exp(x)))
  }
  else{
    stop("Please give a valid target:
         'Gaussian'
         'Logistic'
         ")
  }
  
}
"

MultGauss = function(q){
  return(-dmvn(q, mu = rep(0, length(q)), sigma = V, 
               log = TRUE))
}

grad.MultGauss = function(q){
  return(q%*%solve(V))
}
# %*%solve(V)
# == Product of Logistics ==

prod_logistic = function(q){
  return(sum(-log(exp(q)/(1+exp(q))^2)))
}

grad.prod_logistic = function(q){
  return(-1 + 2*exp(q)/(1+exp(q)))
}


# ==== Kinetic Energy Function ====
# Acts as a proposal in HMC
squared.kinetic = function(rho,m){
  return(sum(0.5*rho^2/m))
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
