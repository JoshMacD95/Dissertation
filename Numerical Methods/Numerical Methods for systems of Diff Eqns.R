
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


Euler = function(q0, rho0, m, timestep, obs.time, grad.U){
  
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
  q.values = matrix(q0, nrow = 1, ncol = d)
  rho.values = matrix(rho0, nrow = 1, ncol = d)
  
  i = 1
  while(i*timestep < obs.time){
    # Adding empty row to fill with new position and momentum values
    q.values = rbind(q.values, rep(0,d))
    rho.values = rbind(rho.values, rep(0,d))
    for (j in 1:d){
    rho.values[i+1,j] = rho.values[i,j] - timestep*grad.U(q.values[i,j])
    q.values[i+1,j] = q.values[i,j] + timestep*rho.values[i,j]/m[j]
    }
    i = i + 1
  }
  return(cbind(q.values, rho.values))
}

# ===================
# ==== Modified Euler
# ===================

Euler.mod = function(q0, rho0, m, timestep, obs.time, grad.U){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensions
  d = length(q0)
  
  q.values = matrix(q0, nrow = 1, ncol = d)
  rho.values = matrix(rho0, nrow = 1, ncol = d)
  i = 1
  while(i*timestep < obs.time){
    # Adding empty row to fill with new position and momentum values
    q.values = rbind(q.values, rep(0, d))
    rho.values = rbind(rho.values, rep(0,d))
    for(j in 1:d){
      rho.values[i+1,j] = rho.values[i,j] - timestep*grad.U(q.values[i,j])
      q.values[i+1,j] = q.values[i,j] + timestep*rho.values[i+1,j]/m[j]
    }
    i = i + 1
  }
  return(cbind(q.values, rho.values))
}

#=====================
# ==== Leapfrog method
#=====================

leapfrog = function(q0, rho0, m, timestep, obs.time, grad.U){
  
  # Making sure mass is defined in each 'dimension'
  if(length(m) == 1){
    m = rep(m, length(rho0))
  }
  if(length(m) != length(rho0)){
    return(print("Error: mass not specified for every direction of motion"))
  }
  
  # d - dimensional problem
  d = length(q0)
  
  q.values = matrix(q0, nrow = 1, ncol = d)
  rho.values = matrix(rho0, nrow = 1, ncol = d)
  
  i = 1
  while(i*timestep < obs.time){
    # Adding empty row to fill with new position and momentum values
    q.values = rbind(q.values, rep(0, d))
    rho.values = rbind(rho.values, rep(0,d))
    for(j in 1:d){
      rho_halfeps = rho.values[i,j] - 0.5*timestep*grad.U(q.values[i,j])
      q.values[i+1,j] = q.values[i,j] + timestep*rho_halfeps/m[j]
      rho.values[i+1,j] = rho_halfeps- 0.5*timestep*grad.U(q.values[i+1,j])
    }
    i = i + 1
  }
  return(cbind(q.values, rho.values))
}

# =====================
# ==== General Function
# =====================

numerical.method = function(q0, rho0, m, timestep, obs.time, grad.U, method, print = FALSE, final = TRUE){
  d = length(q0)
  # Run the choice of numerical method
  result = method(q0, rho0, m, timestep, obs.time, grad.U)
  
  # Store the final position and momentum values
  q.star = tail(result, n = 1)[,1:d]
  rho.star = tail(result, n = 1)[, (d+1):(2*d)]
  
  
  # Following code prints the information about the results, if
  # user asks for it. Reverse calculation will only be carried out if
  # printing is demanded (to remove unneccessary computation)
  if(print == TRUE){
    
    # Feed the final values back into the chosen method, with momentum being reversed.
    reverse = tail(method(q0 = q.star, rho0 = -rho.star, m, timestep, obs.time, grad.U), n = 1)
    
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
    output = c(q.star, rho.star) 
  }
  else{
    output = result
  }
  return(output)
}






