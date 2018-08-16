#'
#' Investigating Optimal Acceptance Rate as Dimension increases
#'
#'
#'


# ==== Preamble ==== 
setwd("~/Dissertation/")
source("HMC Algorithm/Multivariate HMC Algorithm.R")
library(utils)
# ==== Test ====

target = "Prod.Logistic"
method = leapfrog
no.its = 10000
burn.in = 1000

obs.time = 2.5
L.list = 1:20
stepsize.list = obs.time/L.list
d.list = rep(c(500, 1000, 10000), 10)


# Storage
ESS = rep(NA, length(L.list)*length(d.list))
ESS.L = rep(NA, length(L.list)*length(d.list))
AR = rep(NA, length(L.list)*length(d.list))
L = rep(NA, length(L.list)*length(d.list))
Stepsize = rep(NA, length(L.list)*length(d.list))
Dimension = rep(NA, length(L.list)*length(d.list))


# create progress bar
total = length(d.list)*length(L.list) 
pb = txtProgressBar(min = 0, max = total, style = 3)


# Algorithm
for(j in 1:length(d.list)){
  # Offset parameter so values are stored correctly
  k = (j-1)*length(L.list)
  
  for(i in 1:length(L.list)){
    # Start in Stationary Distribution
    theta = 1*((1:d.list[j])%%2 != 0) + 10*((1:d.list[j])%%2 == 0)
    x0 = rlogis(d.list[j], location = 0, scale = 1/theta)
    m = theta
    # Run HMC Algorithm
    HMC = Multivariate.HMC(x0, m, L.list[i], obs.time, target, K = squared.kinetic, method, no.its, burn.in)
    
    # Store Values
    ESS[k + i] = HMC$ESS
    ESS.L[k + i] = HMC$scaled.ESS
    AR[k + i] = HMC$accept.rate
    L[k + i] = L.list[i]
    Stepsize[k + i] = stepsize.list[i]
    Dimension[k + i] = d.list[j]
    setTxtProgressBar(pb, k+i)
  }
}

output.data = data.frame(Dimension, L, Stepsize, ESS, ESS.L, AR)

# == Save Output ==
write.csv(output.data, file = "Output Data/Logistic_HMC_With_Increasing_Dimension3.csv")

