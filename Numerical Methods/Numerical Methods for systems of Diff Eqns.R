
#' Further Topics in MCMC
#' Euler and Leapfrog Numerical Methods
#'
#'
#'
#'

# ==== Hamiltonian Dynamics Equations

# Potential Energy Function (-log(target))
# Assume a Guassian target for now.
U = function(q){
  return(-log(dnorm(0,1)))
}

U_dash = function(q){
  return(q)
}

# Kinetic Energy Function
# Acts as a proposal in HMC
K = function(rho,m){
  return(0.5*rho^2/m)
} 

# Total Energy (This should stay constant by the conservation of energy)
H = function(q, rho, m){
  return(U(q) + K(rho,m))
}



# Hamiltonian Dynamics give rise to a pair of deterministic
# differential equations which typically cannot be solved analytically. 
#Therefore we must use numerical methods.


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


Euler = function(q0, rho0, m, timestep, obs.time){
  
  q.values = c(q0)
  rho.values = c(rho0)
  
  i = 1
  while(i*timestep < obs.time){
    rho.values[i+1] = rho.values[i] - timestep*U_dash(q.values[i])
    q.values[i+1] = q.values[i] + timestep*rho.values[i]/m
    i = i + 1
  }
  return(cbind(q.values, rho.values))
}


# ===================
# ==== Modified Euler
# ===================

Euler.mod = function(q0, rho0, m, timestep, obs.time){
  
  q.values = c(q0)
  rho.values = c(rho0)
  
  i = 1
  while(i*timestep < obs.time){
    rho.values[i+1] = rho.values[i] - timestep*U_dash(q.values[i])
    q.values[i+1] = q.values[i] + timestep*rho.values[i+1]/m
    i = i + 1
  }
  return(cbind(q.values, rho.values))
}



#=====================
# ==== Leapfrog method
#=====================

leapfrog = function(q0, rho0, m, timestep, obs.time){
  
  q.values = c(q0)
  rho.values = c(rho0)
  
  i = 1
  while(i*timestep < obs.time){
    rho_halfeps = rho.values[i] - 0.5*timestep*U_dash(q.values[i])
    q.values[i+1] = q.values[i] + timestep*rho_halfeps/m
    rho.values[i+1] = rho_halfeps- 0.5*timestep*U_dash(q.values[i+1])
    i = i + 1
  }
  return(cbind(q.values, rho.values))
  
}
leapfrog(q0 = 1, rho0 = 0, m =1, timestep = 0.3, obs.time =2*pi)
