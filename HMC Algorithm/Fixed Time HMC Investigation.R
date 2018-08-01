#'
#' Further Topics in MCMC
#' Multivariate HMC Algortihm Investigation (Set Computation Time)
#'
#'
#'




# ==== Preamble ====

library(dplyr)

source("HMC Algorithm/Fixed Time Multivariate HMC Algorithm.R")


# ==== Replicating Beskos ====

d = 1000
x0 = rep(0, d)
m = 1
L = 10
stepsize.list = rep(seq(from = 0.1, to = 1.9, length = 7), 50)
target = "Std.Gaussian"
K = squared.kinetic
comp.time = 30
burn.in = 1000

ESS = c()
Accept.rate = c()
SE = c()
for(i in 1:length(stepsize.list)){
  HMC = HMC.fixedtime(x0, m, L, obs.time = stepsize.list[i]*L, target, K, method = leapfrog, comp.time, burn.in)
  ESS[i] = HMC$ESS
  Accept.rate[i] = HMC$accept.rate
  SE[i] = (mean(HMC$sample[,1]) - 0)^2
}



plot(stepsize.list, Accept.rate)

data = data.frame(cbind(stepsize.list, ESS, Accept.rate))


data = arrange(.data = data, group = stepsize.list)

median.accept.rates = c()

for(i in 1:(length(Accept.rate)/50)){
  median.accept.rates[1:50 + 50*(i-1)] = median(data$Accept.rate[1:50 + 50*(i-1)])
}



data$median.accept.rates = round(median.accept.rates, digits = 3)





boxplot(SE ~ median.accept.rates, data = data,ylim = c(0,0.03),
        ylab = 'Squared Error', xlab = 'Median Acceptance Rate')

boxplot(ESS ~ median.accept.rates, data = data,ylim = c(0,8000),
        ylab = 'ESS', xlab = 'Median Acceptance Rate')



MSE = c()

for(i in 1:(length(SE)/50)){
  MSE[i] = mean(SE[1:50 + 50*(i-1)])
}


median.accept.rates = c()

for(i in 1:(length(Accept.rate)/50)){
  median.accept.rates[i] = median(data$Accept.rate[1:50 + 50*(i-1)])
}



plot(median.accept.rates, MSE, ylim = c(0,0.1))





