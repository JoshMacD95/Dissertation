#'
#' Further Topics in MCMC
#' Effective Sample Size for varying values of stepsize and 
#' Observation Time
#'
#'

# ==== Preamble ====
source("HMC Algorithm/Univariate HMC Algorithm.R")

# ==== Parameters ====

q0 = 0
m = 1
L = 15.7
obs.time = pi/2
U = Gaussian
grad.U = grad.Gaussian
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

# ==== HMC Algorithm ====

test.HMC = Univariate.HMC(q0, m, L, obs.time, 
                          U, grad.U, K, 
                          method = leapfrog, no.its, burn.in)

# ==== Diagnostics ====

# == Stationarity and Auto-Correlation Checking
par(mfrow = c(1,2))

plot(test.HMC$sample, type = 'l') # Check for mixing and convergence

acf(test.HMC$sample) # Check dependence in chain

# == Used to check stationary distribution is correct for normal data

#qqnorm(test.HMC$sample) # If Data is normal, use this to check distribution

#abline(a = 0, b = 1) # The qqplot should follow this line if data is normal

# Plots Histogram of samples and compares with the density of the target distribution      
hist(test.HMC$sample, freq = FALSE)
x = seq(from = -10, to = 10, length = 100)
density.logistic = exp(-Logistic(x))
lines(x, density.logistic, col = 'blue')       

# == Effective Sample Size ==

obs.times = c(pi/2, pi/4, pi/6, pi/8, pi/10)

no.steps = c(1, 2, 4, 8, 16, 32, 64)

ESS.matrix = matrix(nrow = length(obs.times), ncol = length(no.steps), data = 0)

for(i in 1:length(obs.times)){
  for(j in 1:length(no.steps))
  ESS.matrix[i,] = Univariate.HMC(q0, m, L = no.steps[j], obs.time = obs.times[i], 
                          U, grad.U, K, 
                          method, no.its, burn.in, output.type ="ESS")
}

ESS.per.sec.matrix = matrix(nrow = length(obs.times), ncol = length(no.steps), data = 0)

for(i in 1:length(obs.times)){
  ESS.per.sec.matrix[i,] = sapply(X = no.steps, FUN = Univariate.HMC, obs.time = obs.times[i], q0, m, 
                          U, grad.U, K, 
                          method, no.its, burn.in, output.type ="ESS/sec")
}

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
