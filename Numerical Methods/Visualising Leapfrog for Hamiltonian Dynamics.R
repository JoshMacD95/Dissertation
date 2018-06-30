#'
#' Further Topics in MCMC
#' Visualising Leapfrog for Hamiltonian Dynamics
#'

# ==== Preamble ====

source("Numerical Methods/Numerical Methods for systems of Diff Eqns.R")

# ==== Leapfrog ====
# Visualisation is hard for d > 2, so best if d = 1,2
q0 = c(-1.5,-1.55)
rho0 = c(-1, 1)
m = 1
L = 200
obs.time = L*0.02
U = MultGauss
grad.U = grad.MultGauss
K = squared.kinetic

leapfrog.test = leapfrog(q0, rho0, m, L, obs.time, grad.U)

# ==== Path of Position/Momentum Variables ====
par(mfrow = c(1,3))

# Positions
x1.values = seq(-2, 2, length = 1000)
x2.values = seq(-2, 2, length = 1000)

U.dist.values = matrix(NA, ncol = 1000, nrow = 1000)

for(i in 1:1000){
  for(j in 1:1000){
    U.dist.values[i,j] = exp(-U(c(x1.values[i], x2.values[j])))
  }
}

contour(x1.values, x2.values, U.dist.values)
points(leapfrog.test[,1], leapfrog.test[,2], type ='l', pch = 16, col = 'blue')

# Momentum

rho1.values = seq(-2,2, length = 1000)
rho2.values = seq(-2,2, length = 1000)

K.dist.values = matrix(NA, ncol = 1000, nrow = 1000)

for(i in 1:1000){
  for(j in 1:1000){
    K.dist.values[i,j] = exp(-K(c(rho1.values[i], rho2.values[j]), m))
  }
}

contour(rho1.values, rho2.values, K.dist.values)
points(leapfrog.test[,3], leapfrog.test[,4], type = 'b', pch = 16, col ='purple')

# ==== Hamiltonian Values ====

H.values = rep(NA, nrow(leapfrog.test))

for(i in 1:nrow(leapfrog.test)){
  H.values[i] = U(leapfrog.test[i, 1:2]) + K(leapfrog.test[i, 3:4], m)
}

plot(H.values, type = 'b', pch = 16, col = 'palegreen4')
abline(H.values[1], b = 0, lty = 2)

