#'
#'
#'
#' MALA Investigation
#'
#' If L = 1, HMC is reduced to MALA. Hence, we can use the
#' HMC Algorithm to investigate Optimal Acceptance rate of 
#' MALA.

# ==== Preamble ==== 

source("HMC Algorithm/Multivariate HMC Algorithm.R")

# ==== Test ====


m = 1
#target = "Increasing.Scale.Gauss"
target = "Std.Gaussian"
method = leapfrog
no.its = 10000
burn.in = 1000
d = 2
L = 1
sigma = 1
stepsize.list = seq(from = 0.1, to = 2*sigma, by = 0.1)
ESS = rep(NA, length(stepsize.list))
AR = rep(NA, length(stepsize.list))
stepsize = stepsize.list

for(i in 1:length(stepsize.list)){
  x0 = rnorm(d, mean = 0, sd = sigma)
  HMC = Multivariate.HMC(x0, m, L , obs.time = stepsize.list[i], target, K = squared.kinetic, method, no.its, burn.in)
  ESS[i] = HMC$ESS
  AR[i] = HMC$accept.rate
}

dim = rep(2, length(stepsize))
HMC.output = data.frame(dim, stepsize, ESS, AR)

plot(HMC.output$AR, HMC.output$ESS, xlim = c(0,1))

write.csv(HMC.output, file ="MALA_2Dim.csv")
