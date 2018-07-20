#'
#' Further Topics in MCMC
#' Functions for MCMCs (RWM, Independence, Gibbs)
#' 
#'

# ==== Preamble ====

library(mvnfast)

# ==== Targets/Proposals ====

target = function(x, dist){
  if(dist =="1D.Gaussian"){
    return(dnorm(x, mean = 0, sd = 1))
  }
  
  if(dist == "Gaussian"){
    return(dmvn(x, mu = rep(0, length(x)), sigma = diag(1, length(x))))
  }
  
  else if(dist == "1D.t"){
    return(dt(x, df = 4, ncp = 0))
  }
  
  else if(dist == "t"){
    return(dmvt(x, mu = rep(0, length(x)), sigma = diag(1, length(x)), df = 4))
  }
  else if(dist == "Cauchy"){
    return(dcauchy(x, location = 0, scale = 1))
  }
  else{
    stop("Please give a valid target/proposal:
         '1D.Gaussian'
         'Gaussian'
         '1D.t
         't'
         ")
  }
}

proposal = function(x, dist){
  
  if(dist =="1D.Gaussian"){
    return(rnorm(x, mean = 0, sd = 1))
  }
  
  if(dist == "Gaussian"){
    return(rmvn(x, mu = rep(0, length(x)), sigma = diag(1, length(x))))
  }
  
  else if(dist == "1D.t"){
    return(rt(x, df = 4, ncp = 0))
  }
  
  else if(dist == "t"){
    return(rmvt(x, mu = rep(0, length(x)), sigma = diag(1, length(x)), df = 4))
  }
  else if(dist == "Cauchy"){
    return(rcauchy(x, location = 0, scale = 1))
  }
  else{
    stop("Please give a valid target/proposal:
         '1D.Gaussian'
         'Gaussian'
         '1D.t
         't'
         ")
  }
}

output.RWM = function(target, proposal, lambda, ESS, accept.rate, x0){
  # Problem with Target and Kinetic function names
  print(paste(c("Target: ", target), sep = "", collapse = ""))
  print(paste(c("Proposal:", proposal)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size:", ESS)), sep = "", collapse = "")
  print(paste(c("Acceptance Rate:", accept.rate)), sep = "", collapse = "")
  print(paste(c("Initial Conditions:", x0)), sep = "", collapse = "")
  print(paste(c("Number of Iterations:", no.its)), sep = "", collapse = "")
  print(paste(c("Scale Parameter", lambda)), sep = "", collapse = "")
}









