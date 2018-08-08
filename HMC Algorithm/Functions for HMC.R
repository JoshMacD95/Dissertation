#'
#' Further Topics in MCMC
#' Functions for Univariate HMC
#'

#'
#' This file is to be used to store any functions/distributions
#' to be used in Hamiltonian Monte Carlo Algoritms/Numerical Methods.
#'
#' 
#'
#'
# ==== Preamble ====
library(mvnfast)

# function for calculation of exp(x)/(1+exp(x))^2

exp.over.1pexp = function(x){
  y = exp(x)/(1+exp(x))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] = 0
    }
  }
  return(y)
} 

# function for calculation of log(1 + e^x)
log1pexp = function(x){
  y = c()
  for(i in 1:length(x)){
    if(x[i] <= 37){
      y[i] = exp(x[i])
    }
    else if(37<x[i] & x[i]<= 18){
      y[i] = log1p(exp(x[i]))
    }
    else if(18 < x[i] & x[i] <= 33.3){
      y[i] = x[i] + exp(-x[i])
    }
    else{
      y[i] = x[i]
    }
  }
 return(y)
}

# ==== Potential Energy Function ==== 

# NOTE: If a distribution is edited in the U (potential) function
#       please change the corresponding distribution in the grad.U 
#       (Grad Potential) to line up 
U = function(x, target){
  if(target == "1D.Gaussian"){
    return(-dnorm(x, mean = 0, sd = 1, log = TRUE))
  }
  if(target == "Std.Gaussian"){
    return(-dmvn(x, mu = rep(0, length(x)), sigma = diag(1^2, length(x)), log = TRUE))
  }
  else if(target == "General.Gaussian"){
    # Covariance Matrix
    V = diag(1, length(x))
    return(dmvn(x, mu = rep(0, length(x)), sigma = V))
  }
  else if(target == "Logistic"){
    theta = 10
    return(-log(theta) - theta*x + 2*log(1 + exp(theta*x)) )
  }
  #log1pexp(theta*x)
  #else if(target == "Prod.Logistic"){
  #  i = seq(1, length(x))
  #  theta = 1*(i%%2 == 0) + 1*(i%%2 != 0)
  #  return(sum(log(theta) + x*theta - 2*log1pexp(x*theta)))
  #}
  else if(target == "Prod.Logistic"){
    i = seq(1, length(x))
    theta = 10*(i%%2 == 0) + 1*(i%%2 != 0)
    #theta = rep(1, length(x))
    return(sum(-log(theta) - theta*x + 2*log(1 + exp(theta*x))))
  }
  else if(target == "Increasing.Scale.Gauss"){
    return(-dmvn(x, mu = rep(0, length(x)), sigma = diag(seq(from = 0.1, to = length(x)*0.1, length = length(x))), log = TRUE))
  }
  else{
    stop("Please give a valid target:
         '1D.Gaussian'
         'Std.Gaussian'
         'General.Gaussian'
         'Logistic'
         'Prod.Logisitic'
         'Increasing.Scale.Gauss'
         ")
  }
}

#rep(c(1,10), length(x)/2)

grad.U = function(x, target){
  if(target == "1D.Gaussian"){
    return(x)
  }
  else if(target == "Std.Gaussian"){
    return(x/(1)^2)
  }
  else if(target == "General.Gaussian"){
    V = diag(1, length(x))
    return(x%*%solve(V))
  }
  else if(target == "Logistic"){
    theta = 10
    return(-theta + 2*theta*exp(theta*x)/(1+exp(theta*x)))
  }
   
  else if(target == "Prod.Logistic"){
    i = seq(1, length(x))
    theta = 10*(i%%2 == 0) + 1*(i%%2 != 0)
    #theta = rep(1, length(x))
    return(-theta + 2*theta*exp(theta*x)/(1+exp(theta*x)))
  }
  else if(target == "Increasing.Scale.Gauss"){
    return(x%*%solve(diag(seq(from = 0.1, to = length(x)*0.1, length = length(x)))))
  }
  else{
    stop("Please give a valid target:
         '1D.Gaussian'
         'Std.Gaussian'
         'General.Gaussian'
         'Logistic'
         'Prod.Logisitic'
         'Increasing.Scale.Gauss'
         ")
  }
}

# ==== Kinetic Energy Function ====
# Acts as a proposal in HMC
squared.kinetic = function(rho,m){
  return(sum(0.5*rho^2/m))
} 

# ==== Output for HMC Algorithm ====

output = function(target, K, ESS, accept.rate, x0, no.its, L, obs.time){
  # Problem with Target and Kinetic function names
  print(paste(c("Target: ", target), sep = "", collapse = ""))
  print(paste(c("Kinetic Energy Function:", K)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size:", min(ESS))), sep = "", collapse = "")
  print(paste(c("Acceptance Rate:", accept.rate)), sep = "", collapse = "")
  print(paste(c("Initial Conditions:", x0)), sep = "", collapse = "")
  print(paste(c("Number of Iterations:", no.its)), sep = "", collapse = "")
  print(paste(c("Steps", L)), sep = "", collapse = "")
  print(paste(c("Observation Time", obs.time)), sep = "", collapse = "")
}


# Total Energy (This should stay constant by the conservation of energy)
#H = function(q, rho, m){
#  return(U(q) + K(rho,m))
#}

