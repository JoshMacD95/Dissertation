#' 
#' Further Topics in MCMC
#' Multivariate HMC Diagnostics
#'

# ==== Preamble ====

  source("HMC Algorithm/Multivariate HMC Algorithm.R")

# See Radford Neal pd 135 onwards

# ==== What Stepsize? ====
# Depends on Hamiltonian and what makes the hamiltonian unstable
# For Multidimensional problems stability of stepsize will be determined by th
# width of the ditribution in the most constrained direction
# For 1-dimension normal (orthogonal), stepsize < 2*sigma is stable
d = 4
x0 = rnorm(d, 0, 1)
m =1
rho0 = rmvn(1, mu = rep(0,d), sigma = diag(m,d))
L = 1000
stepsize = 0.01
target = "Std.Gaussian"
K = squared.kinetic
no.its = 10000
burn.in = 1000


stepsize.seq = c(0.01, 0.1, 1, 1.5, 1.9)
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


# ==== Premliminary Running ====

d = 4
x0 = rnorm(d, 0, 1)
m = 1
L = 10
#stepsize = 0.155
target = "Std.Gaussian"
K = squared.kinetic
no.its = 10000
burn.in = 1000

Test.MHMC = Multivariate.HMC(x0, m, L, obs.time = 1.4, 
                             target, K, method = leapfrog,
                             no.its, burn.in)
View(Test.MHMC$sample)
par(mfrow = c(1,d))

for(i in 1:d){
  plot(Test.MHMC$sample[,i])
}
for(i in 1:d){
  acf(Test.MHMC$sample[,i], main ='')
}

# ==== Effective Sample Size Vs. Observation Time (varying No. Steps) ==== 

d = 10
x0 = rep(0, d)
m = 1
L.list = seq(2, 10, length = 5)
obs.time.list = seq(from = 0.01, to = pi, length = 10)
target = 'Std.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

scaled.ESS = matrix(NA, nrow = length(obs.time.list), ncol = length(L.list))

for(i in 1:length(obs.time.list)){
  for(j in 1:length(L.list)){
    scaled.ESS[i,j] = Multivariate.HMC(x0, m, L = L.list[j], obs.time = obs.time.list[i],
                                target, K, method, no.its, burn.in)$scaled.ESS
  }
}



plot(obs.time.list, scaled.ESS[,1], col = 1,
     xlab = 'T', ylab = 'ESS/L', ylim = c(min(scaled.ESS), 10000), xlim = c(0,pi), type  = 'l')
for(i in 2:length(L.list)){
  points(obs.time.list, scaled.ESS[,i], col = i+1, type = 'l')
}

legend("topleft", legend = c("2", "4", "6", "8", "10"), title = "Leapfrog Steps",
       col = 1:5, lty = 1)

# ==== Effective Sample Size Vs. Leapfrog Steps (varying Stepsize) ==== 

d = 2
x0 = rep(0, d)
m = 1
L.list = 1:10
stepsize.list = seq(from = 0.01, to = 0.15, length = 5)
target = 'Prod.Logistic'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

ESS = matrix(NA, nrow = length(L.list), ncol = length(stepsize.list ))

for(i in 1:length(L.list)){
  for(j in 1:length(stepsize.list)){
    ESS[i,j] = Multivariate.HMC(x0, m, L = L.list[i], obs.time = stepsize.list[j]*L.list[i],
                                       target, K, method, no.its, burn.in)$ESS
  }
}



plot(L.list, ESS[,1], col = 1,
     xlab = 'L', ylab = 'ESS', type = 'l',
     ylim = c(min(ESS), 10000), xlim = c(1,10))
for(i in 2:length(stepsize.list)){
  points(L.list, ESS[,i], col = i+1, type = 'l')
}

legend("topleft", 
       legend = c(as.character(stepsize.list[1]), as.character(stepsize.list[2]), as.character(stepsize.list[3]), as.character(stepsize.list[4]), as.character(stepsize.list[5])),
       title = "Stepsize", col = 1:5, lty = 1)


# ==== ESS/sec or ESS/L vs. Acceptance rate (Varying Stepsize, Calculating No. leapfrog steps) ====

d = 20
m = 1
k = 0:5
no.runs = length(k)
stepsize.list = 2^(-k)
L.list = 2^k
stepsize.list = seq(from = 0.03125, to = 2, length = 6)
L.list = ceiling(obs.time/stepsize.list)
target = 'Std.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

reps = 10

all.ESS = rep(0,reps*no.runs)
all.ESS.L = rep(0,reps*no.runs)
all.Acceptance.rates = rep(0,reps*no.runs)

for(j in 1:reps){

  for(i in 1:no.runs){
    ESS = c()
    ESS.L = c()
    Acceptance.rates = c()
    x0 = rnorm(d, 0, 1) # Only sample x0 from dist.n if target is trivial gaurantees starting in stationarity
    HMC =  Multivariate.HMC(x0, m, L = L.list[i], obs.time = stepsize.list[i]*L.list[i],
                            target, K, method, no.its, burn.in)
    ESS[i] = HMC$ESS
    ESS.L[i] = HMC$scaled.ESS
    Acceptance.rates[i] = HMC$accept.rate
  }
  all.ESS[(j-1)*no.runs + (1:no.runs)] = ESS
  all.ESS.L[(j-1)*no.runs + (1:no.runs)] = ESS.L
  all.Acceptance.rates[(j-1)*no.runs + (1:no.runs)] = Acceptance.rates
}


par(mfrow=c(1,1))
plot(all.Acceptance.rates, all.ESS.L)
plot(all.Acceptance.rates, all.ESS)

onehundred.dimensions = data.frame(all.ESS, all.ESS.L, all.Acceptance.rates)

# ==== Finding optimal acceptance rate ====

# Goal: To find a optimal acceptance rate. An acceptance rate which
#       results in low dependence in sample, good mixing and convergence
#       and hence results in an Effective sample size, close to no.its-burn.in.

#       The optimal acceptance rate will do all this while keeping 
#       computation to a minimum. To investigate this I will calculate
#       a scaled ESS which takes into account the number of leapfrog steps
#       taken. This gives a good indication of ESS/sec without the confounding 
#       of computational performance.

# The higher the number of leapfrog steps taken, the more accurate the Hamiltonian Dynamics
# This results in a higher acceptance rate as error in the Hamiltonian is reduced
# due to stepsize being reduced as no. steps increases (Constant trajectory).
# However, the more steps that are taken, the more computational power which is needed
# for one proposal and therefore more time is taken. 
# Hence, the trade-off between the added time due to more steps and the effiency of the 
# algorithm must be considered.
x0 = rep(0,2)
m = 1
#obs.time = 1.55
target = 'Std.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000


L.list = 1:10
stepsize.list = seq(from = 0.01, to = 1.99, length = 10)
scaled.ESS = matrix(NA, nrow = length(L.list), ncol = length(stepsize.list))
acceptance.rates = matrix(NA, nrow = length(L.list), ncol = length(stepsize.list))
j = 1  
for(i in 1:length(L.list)){
  for(j in 1:length(stepsize.list)){
  HMC = Multivariate.HMC(x0, m, L.list[i], obs.time = L.list[i]*stepsize.list[j], target, K, method = leapfrog, no.its, burn.in)
  scaled.ESS[i,j] = HMC$scaled.ESS
  acceptance.rates[i,j] = HMC$accept.rate
 }
}

par(mfrow = c(1,1))
plot(acceptance.rates[,1], scaled.ESS[,1], ylab = 'ESS/L', xlab = 'Acceptance Rate', col = 'red',
     xlim = c(0, 1), ylim = c(0, 9000))

for(i in 2:length(stepsize.list)){
  points(acceptance.rates[,i], scaled.ESS[,i], col = i+1)
}

# == Squared Error of Monte Carlo Estimates for varying stepsize ==
d = 10
x0 = rep(0,d)
m = 1
stepsize.list = seq(from = 0.1, to = 2, length = 10)
L.list = 10
target = 'Std.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000

SE = matrix(nrow = length(stepsize.seq), ncol = 120)

accept.matrix = c(nrow = length(stepsize.seq), ncol = 120)
for(j in 1:120){
  for(i in 1:length(stepsize.list)){
    HMC = Multivariate.HMC(x0, m, L.list, obs.time = L.list*stepsize.list[i], target, K, method = leapfrog, no.its, burn.in)
    SE[i] = ((1/length(HMC$sample[,1]))*sum(HMC$sample[,1]) - 0)^2
    accept.matrix[i] = HMC$accept.rate
  }
}



par(mfrow = c(1,1))
plot(accept.matrix, SE, ylab = 'SE', xlab = 'Acceptance Rate', col = 'red',
     xlim = c(0, 1))

for(i in 2:length(stepsize.list)){
  points(accept.matrix[i,], SE[i,], col = i+1)
}

min(SE)

# ==== Increasing Dimension ====

#' How does the optimal acceptance rate increase with dimension
#' Run the same process as above but for many dimensions and 
#' see how the optimal acceptance rate, with respect to ESS/L
#' evolves with dimension


m = 1
obs.time = 1.55
target = 'Std.Gaussian'
K = squared.kinetic
method = leapfrog
no.its = 10000
burn.in = 1000
L.list = 1:10
stepsize.seq = seq(from = 0.01, to = 1.99, length = 10)
d.seq = c(2, 4, 10, 20, 50, 100)

scaled.ESS = matrix(NA, nrow = length(d.seq), ncol = length(L.list))

acceptance.rates = matrix(NA, nrow = length(d.seq), ncol = length(L.list)*10)

for(j in 1:length(d.seq)){
  for(i in 1:length(L.list)){
    HMC = Multivariate.HMC(x0 = rep(0,d.seq[j]) , m, L.list[i], obs.time, target, K, method = leapfrog, no.its, burn.in)
    scaled.ESS[j,i] = HMC$scaled.ESS
    acceptance.rates[j,i] = HMC$accept.rate
  }
}

scaled.ESS2 = scaled.ESS[,1:100]
acceptance.rates2 = acceptance.rates[,1:100]
plot(acceptance.rates2[1,], scaled.ESS2[1,], ylab = 'ESS/L', xlab = 'Acceptance Rate',
     col = 'red', ylim = c(min(scaled.ESS2), max(scaled.ESS2)), xlim = c(min(acceptance.rates2), max(acceptance.rates2)))

points(acceptance.rates[2,], scaled.ESS[2,],
       col = 'purple')


points(acceptance.rates[3,], scaled.ESS[3,],
       col = 'blue')

points(acceptance.rates[4,], scaled.ESS[4,],
       col = 'green')

points(acceptance.rates[5,], scaled.ESS[5,],
       col = 'orange')



aL.seq = seq(1, 10, by = 1)
eps.seq = seq(0.1, 1, length = 10)
ESS = matrix(0, nrow = length(L.seq), ncol = length(eps.seq))
for(i in 1:length(L.seq)){
  for(j in 1:length(eps.seq)){
    ESS[i,j] = min(Multivariate.HMC(x0, m, L.seq[i], obs.time = L.seq[i]*eps.seq[j], target, K, method = leapfrog, no.its, burn.in)$ESS)
  }
}



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

# ==== Investigation of Optimal Acceptance Rate ====

# Dimensions
d.list = d


