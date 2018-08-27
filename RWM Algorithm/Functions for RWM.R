#'
#'
#' Functions for RWM
#'
#'

# ==== Preamble ====
library(mvnfast)

log.density=function(x,target){
  if(target=="Std.Gaussian"){ # standard normal
    d = length(x)
    ld=dmvn(x, mu = rep(0, d), sigma = diag(1, d), log = TRUE)
  }
  else if(target == "Stretched.Gaussian"){
    d = length(x)
    ld=dmvn(x, mu = rep(0, d), sigma = diag(c(1,10), d), log = TRUE)
  }
  else if(target == "Prod.Logistic"){
    i = seq(1, length(x))
    theta = 10*(i%%2 == 0) + 1*(i%%2 != 0)
    #theta = rep(1, length(x))
    return(-sum(-log(theta) - theta*x + 2*log(1 + exp(theta*x))))
  }
  else if (target=="Laplace") { # 0.5 exp(-|x|)
    ld=-abs(x)/2
  }
  else if (target=="t5") {
    ld=dt(x,df=5,log=TRUE)
  }
  else if (target=="Cauchy") {
    ld=dt(x,df=1,log=TRUE)
  }
  else if (target=="bimodal") {
    ld=log(0.2*dnorm(x,mean=-8,sd=1)+0.8*dnorm(x,mean=6,sd=2))
  }
  else {
    stop("Unknown target type (Gaussian/Laplace/t5/Cauchy/bimodal)")
  }
  
}

#
# Create a proposal with scale parameter 1
#

propose=function(proposal, d, V) {
  if (proposal=="Std.Gaussian") {
    z=rmvn(1, mu = rep(0, d), sigma = V)
  }
  else if (proposal=="Laplace") { # 0.5 exp(-|x|)
    z=sample(c(-1,1),1)*rexp(1)
  }
  else if (proposal=="t5") {
    z=rt(1,df=5)
  }
  else if (proposal=="Cauchy") {
    z=rt(1,df=1)
  }
  else {
    stop("Unknown proposal type (Gaussian/Laplace/t5/Cauchy)")
  }
  return(z)
}



# ==== Output for RWM Algorithm ====

output.RWM = function(target, proposal, lambda, ESS, accept.rate, no.its){
  # Problem with Target and Kinetic function names
  print(paste(c("Target: ", target), sep = "", collapse = ""))
  print(paste(c("Proposal:", proposal)), sep = "", collapse = "")
  print(paste(c("Lambda:", lambda)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size:", min(ESS))), sep = "", collapse = "")
  print(paste(c("Acceptance Rate:", accept.rate)), sep = "", collapse = "")
  #print(paste(c("Initial Conditions:", x0)), sep = "", collapse = "")
  print(paste(c("Number of Iterations:", no.its)), sep = "", collapse = "")
}




