#'
#' Further Topics in MCMC
#' Order of Expected Hamiltonian Errors
#'
#'

# ==== Preamble ====

source("HMC Algorithm/Multivariate HMC Algorithm.R")
source("HMC Algorithm/Univariate HMC Algorithm.R")


# ==== Expected Hamiltonian Error is order eps^4 ===

stepsize.list = seq(0.01, 2, by = 0.02)
L = 10
avg.ham.errors = rep(NA, length(stepsize.list))


for(i in 1:length(stepsize.list)){
  HMC = Univariate.HMC(rnorm(1,0,1), m = 1, L, obs.time = stepsize.list[i]*L, target = "1D.Gaussian", K = squared.kinetic, no.its = 10000, burn.in = 1000)
  avg.ham.errors[i] = mean(HMC$delta.H)
}
rnorm(0,1)

par(mfrow = c(1,2))
plot(stepsize.list, avg.ham.errors)
plot(stepsize.list, avg.ham.errors^(1/4))
eps = seq(0.001, 1, by = 0.001)
y = eps 
lines(eps, y, col = 'green', type ='l')

# ==== Expected Hamiltonian Error grows linearly with dimension ====

stepsize.list = seq(0.01, 1, by = 0.02)

plot(n, eps.n)










