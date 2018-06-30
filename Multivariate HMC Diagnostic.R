#' 
#' Further Topics in MCMC
#' Multivariate HMC Diagnostics
#'

# ==== Preamble ====

source("HMC Algorithm/Multivariate HMC Algorithm.R")

# ==== Algorithm ====

q0 = c(0,0)
m = 1
L = 25
obs.time = 2*pi
U = MultGauss
grad.U = grad.MultGauss
K = squared.kinetic
no.its = 10000
burn.in = 1000
d = length(q0)

Test.MHMC = Multivariate.HMC(q0 = c(0,0), m, L, obs.time, 
                             U, grad.U, K, method = leapfrog,
                             no.its, burn.in, output.type = "all")

par(mfrow = c(1,1))
plot(Test.MHMC$hamiltonian, lty=2, type = 'l')

# ==== Mixing, Convergence and Chain Dependence ====

par(mfrow= c(2,d))

for(i in 1:d){
  plot(Test.MHMC$sample[,i], type = 'l')
  acf(Test.MHMC$sample[,i])
}

# ==== Contour Plotting ====

x_1 = seq(from = -5, to = 5, length = 1000)
x_2 = seq(from = -5, to = 5, length = 1000)
mvn.values = matrix(data = NA, nrow = 1000, ncol = 1000)

for(i in 1:1000){
  for(j in 1:1000){
    mvn.values[i,j] = dmvn(X = c(x_1[i], x_2[j]), mu = rep(0, 2), sigma = matrix(data = c(1, 0, 0, 1), ncol = 2, byrow =T))
  }
}

par(mfrow=c(1,1))
contour(x_1, x_2, mvn.values)
points(Test.MHMC$sample[1:1000,1], Test.MHMC$sample[1:1000,2], pch = 4)
