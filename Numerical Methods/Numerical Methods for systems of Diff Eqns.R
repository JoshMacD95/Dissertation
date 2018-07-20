
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
x_t = function(x0, m, t){
  return(x0*cos(t/sqrt(m)))
}

rho_t = function(x0, m, t){
  return(-(x0/sqrt(m))*sin(t/sqrt(m)))
}


# ==========
# ==== Euler
# ==========


Euler = function(x0, rho0, m, L, obs.time, target, K){

  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensions
  d = length(x0)
  
  # Creating Empty Matrix for Values
  x.values = matrix(0, nrow = L + 1, ncol = d)
  rho.values = matrix(0, nrow = L + 1, ncol = d)
  hamiltonian = c()
  # Intial Values
  x.values[1,] = x0
  rho.values[1,] = rho0
  hamiltonian[1] = U(x0, target) + K(rho0)
  for(i in 1:L){
    rho.values[i+1,] = rho.values[i,] - (obs.time/L)*grad.U(x = x.values[i,], target = target)
    x.values[i+1,] = x.values[i,] + (obs.time/L)*rho.values[i,]/m
    hamiltonian[i+1] = U(x.values[i+1,], target) + K(rho.values[i+1,], m)
  }
  return(list(x.values = x.values, rho.values = rho.values, hamiltonian = hamiltonian))
}

# ===================
# ==== Modified Euler
# ===================

Euler.mod = function(x0, rho0, m, L, obs.time, target, K){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensions
  d = length(x0)
  
  # Creating empty matrix for values
  x.values = matrix(0, nrow = L + 1, ncol = d)
  rho.values = matrix(0, nrow = L + 1, ncol = d)
  hamiltonian = c()
  # Initial Values
  x.values[1,] = x0
  rho.values[1,] = rho0
  hamiltonian[1] = U(x0, target) + K(rho0, m)
  
  for(i in 1:L){
    rho.values[i+1,] = rho.values[i,] - (obs.time/L)*grad.U(x.values[i,], target)
    x.values[i+1,] = x.values[i,] + (obs.time/L)*rho.values[i+1,]/m
    hamiltonian[i+1] = U(x.values[i+1,], target) + K(rho.values[i+1,], m)
  }
  return(list(x.values = x.values, rho.values = rho.values, hamiltonian = hamiltonian))
}

#=====================
# ==== Leapfrog method
#=====================

# Give number os steps instead o epsilon
leapfrog = function(x0, rho0, m, L, obs.time, target, K){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensional problem
  d = length(x0)
  
  # Creating empty matrix for values
  x.values = matrix(0, nrow = L + 1, ncol = d)
  rho.values = matrix(0, nrow = L + 1, ncol = d)
  hamiltonian = c()
  # Initial Values
  x.values[1,] = x0
  rho.values[1,] = rho0
  hamiltonian[1] = U(x0, target) + K(rho0, m)
  for(i in 1:L){
    rho_halfeps = rho.values[i,] - 0.5*(obs.time/L)*grad.U(x.values[i,], target)
    x.values[i+1,] = x.values[i,] + (obs.time/L)*rho_halfeps/m
    rho.values[i+1,] = rho_halfeps - 0.5*(obs.time/L)*grad.U(x.values[i+1,], target)
    hamiltonian[i+1] = U(x.values[i+1,], target) + K(rho.values[i+1,], m)
  }
  return(list(x.values = x.values, rho.values = rho.values, hamiltonian = hamiltonian))
}






# =====================
# ==== General Function
# =====================

numerical.method = function(x0, rho0, m, L, obs.time, target, K, method, print = FALSE, final = TRUE){
  d = length(x0)
  # Run the choice of numerical method
  result = method(x0, rho0, m, L, obs.time, target, K)
                 
  
  # Store the final position and momentum values
  x.star = result$x.values[L+1,]
  rho.star = result$rho.values[L+1,]
  
  
  # Following code prints the information about the results, if
  # user asks for it. Reverse calculation will only be carried out if
  # printing is demanded (to remove unneccessary computation)
  if(print == TRUE){
    
    # Feed the final values back into the chosen method, with momentum being reversed.
    reverse = tail(method(x0 = x.star, rho0 = -rho.star, m, L, obs.time, grad.U), n = 1)
    
    #print(paste(c("Method:", method.name)), sep = "", collapse = "")
    print(paste(c("Starting Point:", x0, rho0)), sep = "", collapse = "")
    print(paste(c("Result:", x.star, rho.star)), sep = "", collapse = "")
    print(paste(c("Reverse:", reverse[,1:d], reverse[, (1+d):(2*d)])), sep = "", collapse = "")
  }
  if(final == TRUE){
    output = cbind(x.star, rho.star) 
  }
  else{
    output = result
  }
  return(output)
}





