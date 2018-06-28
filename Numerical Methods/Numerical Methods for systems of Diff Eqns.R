
#' Further Topics in MCMC
#' Numerical Methods for Systems of Differential Equations
#'
#'
#'
#'

# Hamiltonian Dynamics give rise to a pair of deterministic
# differential equations which typically cannot be solved analytically. 
#Therefore we must use numerical methods.


# ==== Preamble ====

source('HMC Algorithm/Functions for HMC.R')




# Analytical Solutions to Hamiltonian Dynamics if possible.
q_t = function(q0, m, t){
  return(q0*cos(t/sqrt(m)))
}

rho_t = function(q0, m, t){
  return(-(q0/sqrt(m))*sin(t/sqrt(m)))
}


# ==========
# ==== Euler
# ==========


Euler = function(q0, rho0, m, L, obs.time, grad.U){

  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensions
  d = length(q0)
  
  #
  q.values = matrix(q0, nrow = L + 1, ncol = d)
  rho.values = matrix(rho0, nrow = L + 1, ncol = d)
  
  for(i in 1:L){
    rho.values[i+1,] = rho.values[i,] - (obs.time/L)*grad.U(q.values[i,])
    q.values[i+1,] = q.values[i,] + (obs.time/L)*rho.values[i,]/m
  }
  return(cbind(q.values, rho.values))
}

Euler(q0 = c(0,0), rho0 = c(1,1), m =1, L = 2, obs.time = 2*pi, grad.U = grad.Gaussian)

# ===================
# ==== Modified Euler
# ===================

Euler.mod = function(q0, rho0, m, L, obs.time, grad.U){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensions
  d = length(q0)
  
  q.values = matrix(q0, nrow = L + 1, ncol = d)
  rho.values = matrix(rho0, nrow = L + 1, ncol = d)

  for(i in 1:L){
    rho.values[i+1,] = rho.values[i,] - (obs.time/L)*grad.U(q.values[i,])
    q.values[i+1,] = q.values[i,] + (obs.time/L)*rho.values[i+1,]/m
  }
  return(cbind(q.values, rho.values))
}

#=====================
# ==== Leapfrog method
#=====================

# Give number os steps instead o epsilon
leapfrog = function(q0, rho0, m, L, obs.time, grad.U){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensional problem
  d = length(q0)
  
  q.values = matrix(q0, nrow = L + 1, ncol = d)
  rho.values = matrix(rho0, nrow = L + 1, ncol = d)

  for(i in 1:L){
    rho_halfeps = rho.values[i,] - 0.5*(obs.time/L)*grad.U(q.values[i,])
    q.values[i+1,] = q.values[i,] + (obs.time/L)*rho_halfeps/m
    rho.values[i+1,] = rho_halfeps - 0.5*(obs.time/L)*grad.U(q.values[i+1,])
  }
  return(cbind(q.values, rho.values))
}






# =====================
# ==== General Function
# =====================

numerical.method = function(q0, rho0, m, L, obs.time, grad.U, method, print = FALSE, final = TRUE){
  d = length(q0)
  # Run the choice of numerical method
  result = method(q0, rho0, m, L, obs.time, grad.U)
  
  # Store the final position and momentum values
  q.star = result[L+1,1:d]
  rho.star = result[L+1,(d+1):(2*d)]
  
  
  # Following code prints the information about the results, if
  # user asks for it. Reverse calculation will only be carried out if
  # printing is demanded (to remove unneccessary computation)
  if(print == TRUE){
    
    # Feed the final values back into the chosen method, with momentum being reversed.
    reverse = tail(method(q0 = q.star, rho0 = -rho.star, m, L, obs.time, grad.U), n = 1)
    
    #if(method == Euler){
    #  method.name = 'Euler'
    #}
    #if(method == Euler.mod){
    #  method.name = 'Modified Euler'
    #}
    #if(method == leapfrog){
    #  method.name = 'Leapfrog'
    #}
    #else{
    #  method.name = 'Not recognised'
    #}
    
    #print(paste(c("Method:", method.name)), sep = "", collapse = "")
    print(paste(c("Starting Point:", q0, rho0)), sep = "", collapse = "")
    print(paste(c("Result:", q.star, rho.star)), sep = "", collapse = "")
    print(paste(c("Reverse:", reverse[,1:d], reverse[, (1+d):(2*d)])), sep = "", collapse = "")
  }
  if(final == TRUE){
    output = cbind(q.star, rho.star) 
  }
  else{
    output = result
  }
  return(output)
}





