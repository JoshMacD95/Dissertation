#'
#' Further Topics in MCMC
#' Functions for MCMCs (RWM, Independence, Gibbs)
#' 
#'

# ==== Preamble ====

library(mvnfast)
library(coda)

# ==== Targets/Proposals ====

density.fn = function(x, target, mu, V, df){
  if(target == "Gaussian"){
    return(dmvn(x, mu, sigma = V))
  }
  
  else if(target == "tdist"){
    return(dmvt(x, mu, sigma = V, df))
  }
  else{
    stop("Please give a valid target/proposal:
         'Gaussian'
         'tdist'
         ")
  }
}







