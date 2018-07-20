#'
#' Further Topics in MCMC
#' Functions for Univariate HMC
#'

# ==== Preamble ====
library(mvnfast)

# ==== Potential Energy Function ==== 

U = function(x, target){
  if(target == "1D.Gaussian"){
    return(-dnorm(x, mean = 0, sd = 1, log = TRUE))
  }
  if(target == "Std.Gaussian"){
    return(-dmvn(x, mu = rep(0, length(x)), sigma = diag(rep(c(1,10), length(x)/2)), 
                 log = TRUE))
  }
  else if(target == "Logistic"){
    return(sum(-log(exp(x)/(1+exp(x))^2)))
  }
  else if(target == "Prod.Logistic"){
    return(sum(-log(exp(x)/(1+exp(x))^2)))
  }
  else{
    stop("Please give a valid target:
         '1D.Gaussian'
         'Std.Gaussian'
         'Logistic'
         'Prod.Logisitic'
         ")
  }
}


grad.U = function(x, target){
  if(target == "1D.Gaussian"){
    return(x)
  }
  if(target == "Std.Gaussian"){
    return(x%*%solve(diag(1,length(x))))
  }
  else if(target == "Logistic"){
    return(-1 + 2*exp(x)/(1+exp(x)))
  }
  else if(target == "Prod.Logistic"){
    return(sum(-1 + 2*exp(x)/(1+exp(x))))
  } 
  else{
    stop("Please give a valid target:
         '1D.Gaussain'
         'Std.Gaussian'
         'Logistic'
         'Prod.Logistic'
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
