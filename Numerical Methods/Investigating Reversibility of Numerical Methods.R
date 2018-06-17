#'
#' Further Topics in MCMC
#' Investingating Reversibility of Numerical Methods
#'

# Given that starting values q0 and rho0 give result q* and rho*, 
# then a starting position of q* and -rho* should give a result of 
# q0, -rho0.

# This code investigates this property of 'reversibility'

# ==== Preamble ====

source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')

# ==== Reversibility ====

# Choose Initial position and momentum variables

q0 = rep(0,50)
rho0 = rep(1,50)
m = 1
timestep = 0.3
obs.time = 2*pi
grad.U = grad.Gaussian

a = numerical.method(q0, rho0, m, timestep, obs.time, grad.U, method=leapfrog, print = FALSE, final = TRUE)


