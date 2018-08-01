#'
#' Further Topics in MCMC
#' Increasing SD Distribution
#'


#'
#'Radford-Neal Constructed a Standard Gaussian Distribution with dimensions of loosing constraint
#'i.e the first dimension has std. dev. 0.01 and the 100th dimension has std. dev. 1
#'
#'The choice of stepsize (epsilon) is constrained by the most constrained direction (the 1st)
#'
#' We still need to move enough in the least constrained direction to propose a move far enough away to reduce
#' auto-correlation. Hence many leapfrogs are needed due to the epsilon being small. 
#'


# ==== Preamble ==== 

source("HMC Algorithm/Multivariate HMC Algorithm.R")

# ==== Test ====

m = 1
#target = "Increasing.Scale.Gauss"
target = "Std.Gaussian"
method = leapfrog
no.its = 10000
burn.in = 1000
d = 50
min.sigma = 0.25

k = 1:5
stepsize.list = 0.425*(2^(-k))
#stepsize.list = seq(2*min.sigma/k, 2*min.sigma, length = k)
L.list = 2^k
#L.list = ceiling(obs.time/stepsize.list)
ESS = rep(0, length(k))
ESS.L = rep(0, length(k))
AR = rep(0, length(k))

for(i in 1:length(k)){
  x0 = rmvn(1, mu = rep(0, d), sigma = diag(min.sigma^2, d))
  #x0 = rmvn(1, mu = rep(0, d), sigma = diag(seq(min.sigma, min.sigma*d, length = d)))
  HMC = Multivariate.HMC(x0, m, L.list[i], obs.time = L.list[i]*stepsize.list[i], target, K = squared.kinetic, method, no.its, burn.in)
  ESS[i] = HMC$ESS
  ESS.L[i] = HMC$scaled.ESS
  AR[i] = HMC$accept.rate
}

dim = rep(d, length(k))

output.data = data.frame(dim, ESS, ESS.L, AR, L.list, stepsize.list)

HMC = Multivariate.HMC(x0, m, L = 15, obs.time = 0.425, target, K = squared.kinetic, method, no.its, burn.in)









