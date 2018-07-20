#'
#' Further Topics in MCMC
#' Univariate HMC Diagnostics
#'
#'

# ==== Preamble ====
library(sna)

source("HMC Algorithm/Univariate HMC Algorithm.R")
# ==== Parameters ====

# ==== What Stepsize? ====
# Depends on Hamiltonian and what makes the hamiltonian unstable
# For Multidimensional problems stability of stepsize will be determined by th
# width of the ditribution in the most constrained direction
# For 1-dimension normal (orthogonal), stepsize < 2*sigma is stable

x0 = 0
m = 1
rho0 = rnorm(1, mean  = 0, sd = sqrt(m))
L = 100
stepsize = 0.01
target = '1D.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000


stepsize.seq = c(0.01, 0.1, 0.5, 1,2)
h = matrix(0, nrow = length(stepsize.seq), ncol = L + 1)
for(i in 1:length(stepsize.seq)){
  h[i,] = leapfrog(x0, rho0, m, L, obs.time = L*stepsize.seq[i], target, K)$hamiltonian
}
par(mfrow =c(1,1))
plot(x = seq(1, length(h[1,])), y = rep(h[1,1], length(h[1,])), 
     type = 'l', col = 'green', lty = 2, 
     ylim= c(min(h),max(h)),
     xlab = "Step Number", ylab = "Hamiltonian"
)
for(i in 1:length(stepsize.seq)){
  lines(h[i,], type='l', col = i)
}

# ==== Preliminary Runs ====

# Goal: find a good observation time which gives good acf and mixing 
#       by varying the number of steps
x0 = 0
m = 1
L = 10
stepsize = 1.55
target = '1D.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000


test.HMC = Univariate.HMC(x0, m, L, obs.time = L*stepsize, 
                          target, K, 
                          method = leapfrog, no.its, burn.in)


par(mfrow = c(1,2))

plot(test.HMC$sample) # Checking for mixing and convergence to stationary distribution
acf(test.HMC$sample) # Checking for dependency in samples

# ==== Finding optimal acceptance rate ====

# Goal: To find a optimal acceptance rate. An acceptance rate which
#       results in low dependence in sample, good mixing and convergence
#       and hence results in an Effective sample size, close to no.its-burn.in.

#       The optimal acceptance rate will do all this while keeping 
#       computation to a minimum. To investigate this I will calculate
#       a scaled ESS which takes into account the number of leapfrog steps
#       taken. This gives a good indication of ESS/sec without the confounding 
#       of computational performance.

x0 = 0
m = 1
obs.time = 1.55
target = '1D.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

L.list = seq(from = 1, to = 50, by = 1)

stepsize.list = obs.time/L.list

scaled.ESS = matrix(nrow = 10, ncol = length(L))
acceptance.rates = 
for(i in 1:length(L.list)){
    HMC = Univariate.HMC(x0, m, L.list[i], obs.time, target, K, method = leapfrog, no.its, burn.in)
    scaled.ESS[i] = HMC$scaled.ESS
    acceptance.rates[i] = HMC$accept.rate
}

par(mfrow = c(1,1))
plot(acceptance.rates[1:20], scaled.ESS[1:20], ylab = 'ESS/L', xlab = 'Acceptance Rate', col = 'red')


# ==== Diagnostics ====

# == Stationarity and Auto-Correlation Checking
par(mfrow = c(1,2))

plot(test.HMC$sample, type = 'l') # Check for mixing and convergence

acf(test.HMC$sample) # Check dependence in chain

# == Used to check stationary distribution is correct for normal data

#qqnorm(test.HMC$sample) # If Data is normal, use this to check distribution

#abline(a = 0, b = 1) # The qqplot should follow this line if data is normal

# Plots Histogram of samples and compares with the density of the target distribution      
par(mfrow = c(1,1))
hist(test.HMC$sample, freq = FALSE)
x.seq = seq(from = -10, to = 10, length = 100)
density = exp(-U(t(x.seq), target))
lines(x, density, col = 'blue')       

# == Effective Sample Size ==

obs.times = c(pi/2, pi/3, pi/4, pi/6, pi/8, pi/10)
rownames = c("pi/2", "pi/3", "pi/4", "pi/6", "pi/8", "pi/10")

no.steps = c(1, 2, 4, 8, 16, 32, 64)
colnames = as.character(no.steps)

ESS.matrix = matrix(nrow = length(obs.times), ncol = length(no.steps), data = 0)

for(i in 1:length(obs.times)){
  for(j in 1:length(no.steps)){
  ESS.matrix[i,j] = Univariate.HMC(x0, m, L = no.steps[j], obs.time = obs.times[i], 
                          target, K, method, no.its, burn.in, 
                          output.type ="ESS")
  }
}

ESS.per.sec.matrix = matrix(nrow = length(obs.times), ncol = length(no.steps), data = 0)

for(i in 1:length(obs.times)){
  for(j in 1:length(no.steps)){
  ESS.per.sec.matrix[i,j] = Univariate.HMC(x0, m, L = no.steps[j], obs.time = obs.times[i],
                                          target, K, method, no.its, burn.in, 
                                          output.type ="ESS/sec")
  }
}

# ==== Effective Sample Size Vs. No. leapsteps (varying stepsize) ==== 

x0 = 0
m = 1
L.list = 1:10
stepsize.list = seq(from = 0.01, to = 1.99, length = 10)
target = '1D.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

scaled.ESS = matrix(NA, nrow = length(L.list), ncol = length(stepsize.list))

for(i in 1:length(L.list)){
  for(j in 1:length(stepsize.list)){
    scaled.ESS[i,j] = Univariate.HMC(x0, m, L = L.list[j], obs.time = stepsize.list[j]*L.list[i],
                                             target, K, method, no.its, burn.in)$ESS
  }
}



plot(L.list, scaled.ESS[,1], col = 1,
     xlab = 'L', ylab = 'ESS')
for(i in 2:length(stepsize.list)){
  plot(L.list, scaled.ESS[,i], col = i+1, type = 'l')
}


plot(no.steps, ESS.per.sec.matrix[1,], 
     type = 'l', col = 1,
     xlab = 'Number of Leapfrog Steps', ylab = 'Effective Sample Size per Second',
     ylim = c(0,max(ESS.per.sec.matrix)))
for(i in 2:length(obs.times)){
  lines(no.steps, ESS.per.sec.matrix[i,], type = 'l', lty = 2, col = i)
}

ESS.matrix = data.frame(ESS.matrix)

colnames(ESS.matrix) = colnames
rownames(ESS.matrix) = rownames

matrixplot(ESS.matrix)
image(ESS.matrix)

# Information about effective sample size
# https://golem.ph.utexas.edu/category/2014/12/effective_sample_size.html


# If obs.time is an even multiple of pi, markov chain will have high, POSITIVE, correlation as 
# proposed values will be close to the current value by hamiltonian dynamics (Full cycles)

# If Obs.time is an odd multiple of pi, markov chain will have high correlation with previous values
# however this correlation alternates between positive and negative 

#  Chosing an obs.time which is between and odd and even multiple of pi will fix these correlation problems
#  giving samples with less correlaation (positive or negative)

#'
#'
#'
#'
#'
#'
#'
# Obs.time too high starts sticking (i.e 100*pi)

# Interesting Result for obs.time = 25*pi
# Explores much faster. Periods of higher variance though
# Seems it occurs for every multiple of 5*pi?

# 15*pi mixes incredibly well, but has moments of sticking.
# Increasing variance of chain

# Problem is Variance parameter, m = 0.5 seems like a good choice. Why?

a = matrix(NA, nrow = 50, ncol = 5)


b = matrix(0 , nrow = 4, ncol = 5)

a[1:4, ] = b