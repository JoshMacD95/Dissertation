#'
#' Further Topics in MCMC
#' Visualising Leapfrog for Hamiltonian Dynamics
#'

# ==== Preamble ====

source("Numerical Methods/Numerical Methods for systems of Diff Eqns.R")

# ==== Contour Plot ====
#target = "Prod.Logistic"
target = "Std.Gaussian"
par(mfrow = c(1,1))

# Positions
x1.values = seq(-2.5, 2.5, length = 1000)
x2.values = seq(-2.5, 2.5, length = 1000)

U.dist.values = matrix(NA, ncol = 1000, nrow = 1000)

for(i in 1:1000){
  for(j in 1:1000){
    U.dist.values[i,j] = exp(-U(c(x1.values[i], x2.values[j]), target))
  }
}


# ==== Leapfrog ====
# Visualisation is hard for d > 2, so best if d = 1,2

x0 = rnorm(2, mean = 0, sd = 1)
m = c(1,1)
rho0 = rmvn(1, mu = rep(0, length(x0)), sigma = diag(m))
L = 10
obs.time = 1.4
K = squared.kinetic

leapfrog.test = leapfrog(x0, rho0, m, L, obs.time, target, K)

# ==== Path of Position/Momentum Variables ====

contour(x1.values, x2.values, U.dist.values)
points(leapfrog.test$x.values[,1], leapfrog.test$x.values[,2], type ='l', pch = 16, col = 'blue')


# ==== Hamiltonian Values ====

H.values = leapfrog.test$hamiltonian

plot(H.values, type = 'b', pch = 16, col = 'palegreen4')

abline(H.values[1], b = 0, lty = 2)

accept.prob = min(c(1, exp(-H.values[length(H.values)]+H.values[1])))
# ==== Path of Momentum variables ====
"
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

"