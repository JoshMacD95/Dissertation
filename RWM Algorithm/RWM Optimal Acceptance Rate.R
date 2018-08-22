#'
#' Investigating Optimal Acceptance Rate as Dimension increases
#'
#'
#'


# ==== Preamble ==== 
source("RWM Algorithm/Multivariate RWM Algorithm .R")
library(utils)
# ==== Test ====

target = "Std.Gaussian"
proposal = "Std.Gaussian"
no.its = 10000
burn.in = 1000

lambda.list = seq(0.001, 3, length = 50)
d.list = rep(12, 10)


# Storage
ESS = rep(NA, length(lambda.list)*length(d.list))
ESS.sec = rep(NA, length(lambda.list)*length(d.list))
AR = rep(NA, length(lambda.list)*length(d.list))
lambda = rep(NA, length(lambda.list)*length(d.list))
Dimension = rep(NA, length(lambda.list)*length(d.list))


# create progress bar
total = length(d.list)*length(lambda.list) 
pb = txtProgressBar(min = 0, max = total, style = 3)


# Algorithm
for(j in 1:length(d.list)){
  # Offset parameter so values are stored correctly
  k = (j-1)*length(lambda.list)
  
  for(i in 1:length(lambda.list)){
    # Start in Stationary Distribution
    #theta = 1*((1:d.list[j])%%2 != 0) + 10*((1:d.list[j])%%2 == 0)
    #x0 = rlogis(d.list[j], location = 0, scale = 1/theta)
    x0 = rnorm(d.list[j])
    # Run HMC Algorithm
    RWM = Multivariate.RWM(no.its, lambda = lambda.list[i], target, proposal, x0, burn.in)
    
    # Store Values
    ESS[k + i] = RWM$ESS
    ESS.sec[k + i] = RWM$scaled.ESS
    AR[k + i] = RWM$accept.rate
    lambda[i] = lambda.list[i]
    Dimension[k + i] = d.list[j]
    setTxtProgressBar(pb, k+i)
  }
}

output.data = data.frame(Dimension, lambda, ESS, ESS.sec, AR)

plot(output.data$AR, output.data$ESS.sec,ylim = c(0,400))
abline(v = 0.234, col = 'red', lty = 2)
# == Save Output ==
write.csv(output.data, file = "Output Data/Optimal_RWMDimension12_2.csv")
