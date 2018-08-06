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


#target = "Increasing.Scale.Gauss"

m = 1
target = "Prod.Logistic"
method = leapfrog
no.its = 1000
burn.in = 10
min.sigma = 1

# How many different combinations per dimension?
d = 2
m = 1
theta = 0.01*((1:d)%%2 != 0) + 0.01*((1:d)%%2 == 0) 
x0 = c(0,0)
HMC = Multivariate.HMC(x0, m, L = 100, obs.time = 100*0.1, target, K = squared.kinetic, method, no.its, burn.in)

plot(HMC$sample[,2])
x0 = rlogis(d, location = 0, scale = 1/theta)

# Integration time, chosen to be adequate but not optimal
obs.time = 1

# No. Leapfrog steps to iterate over
L.list = c(6, 7, 8, 9, 10, 13, 27)
#L.list = ceiling(obs.time/stepsize.list)
#stepsize.list = seq(2*min.sigma/k, 2*min.sigma, length.out = k)

d.list = c(2, 5, 10, 50, 100)

# Stepsize calculated based on Integration time and No. Leapfrog steps
stepsize.mat = matrix(0, nrow = length(L.list), ncol = length(d.list))

for(i in 1:length(L.list)){
  for(j in 1:length(d.list)){
    stepsize.mat[i,j] = obs.time/(L.list[i]*d.list[j]^(1/4))
  }
}

# Storage
ESS = rep(NA, length(L.list)*length(d.list))
ESS.L = rep(NA, length(L.list)*length(d.list))
AR = rep(NA, length(L.list)*length(d.list))
L = rep(NA, length(L.list)*length(d.list))
Stepsize = rep(NA, length(L.list)*length(d.list))
Dimension = rep(NA, length(L.list)*length(d.list))


# Algorithm
for(j in 1:length(d.list)){
  # Offset parameter so values are stored correctly
  k = (j-1)*length(L.list)
  
  for(i in 1:length(L.list)){
    # Start in Stationary Distribution
    x0 = rmvn(1, mu = rep(0, d.list[j]), sigma = diag(min.sigma^2, d.list[j]))
    #x0 = rmvn(1, mu = rep(0, d), sigma = diag(seq(min.sigma, min.sigma*d, length = d)))
    
    # Run HMC Algorithm
    HMC = Multivariate.HMC(x0, m, L.list[i], obs.time = L.list[i]*stepsize.mat[i,j], target, K = squared.kinetic, method, no.its, burn.in)
    
    # Store Values
    ESS[k + i] = HMC$ESS
    ESS.L[k + i] = HMC$scaled.ESS
    AR[k + i] = HMC$accept.rate
    L[k + i] = L.list[i]
    Stepsize[k + i] = stepsize.mat[i,j]
    Dimension[k + i] = d.list[j]
  }
}

output.data = data.frame(Dimension, L, Stepsize, ESS, ESS.L, AR)

plot(output.data$AR, output.data$ESS.L, col = 1, ylim = c(0, 7000), xlim = c(0,1))
for(i in 2:length(Dimension)){
  points(output.data$AR[output.data$Dimension == Dimension[i]], output.data$ESS.L[output.data$Dimension == Dimension[i]], col = i, type = 'l')
}

# == Save Output ==
write.csv(output.data, file = "Dimensionality of HMC")

# == Load Data ==
output = read.csv("Dimensionality of HMC")[,-1]

Dimension = c(2, 5, 10, 20, 50, 100)

plot(output$AR[output$Dimension== 2], output$ESS.L[output$Dimension ==2], col = 1, ylim = c(0, max(output$ESS.L)), xlim = c(0,1))
for(i in 2:length(Dimension)){
  points(output$AR[output$Dimension == Dimension[i]], output$ESS.L[output$Dimension == Dimension[i]], col = i)
}



# Keep observation time below a value which is optimal given ample 
# leapfrog steps.
m = 1
#target = "Increasing.Scale.Gauss"
target = "Std.Gaussian"
method = leapfrog
no.its = 10000
burn.in = 1000
d = 2
sigma = 1

max.obs.time = 1.6
min.obs.time = 0.8
max.stepsize = 2*sigma

stepsize.list = seq(from = 0.1, to = max.stepsize, by = 0.01)
L.list = 1:10

ESS = rep(NA, length(stepsize.list)*length(L.list))
ESS.L = rep(NA, length(stepsize.list)*length(L.list))
AR = rep(NA, length(stepsize.list)*length(L.list))
L = rep(NA, length(stepsize.list)*length(L.list))
stepsize = rep(NA, length(stepsize.list)*length(L.list))
k = 1

for(i in 1:length(stepsize.list)){
  for(j in 1:length(L.list)){
    obs.time = stepsize.list[i]*L.list[j]
    if(min.obs.time < obs.time & obs.time < max.obs.time){
      x0 = rnorm(d, mean = 0, sd = sigma)
      HMC = Multivariate.HMC(x0, m, L = L.list[j], obs.time = obs.time, target, K = squared.kinetic, method, no.its, burn.in)
      ESS[k] = HMC$ESS
      ESS.L[k] = HMC$scaled.ESS
      AR[k] = HMC$accept.rate
    }
    L[k] = L.list[j]
    stepsize[k] = stepsize.list[i]
    k = k + 1
  }
}

dim = rep(d, length(stepsize.list)*length(L.list))

HMC.output = data.frame(dim, L, stepsize, ESS, ESS.L, AR)

write.csv(HMC.output, file ="2_Dimensions_HMC")
table = read.csv("2_Dimensions_HMC")

plot(HMC.output$AR, HMC.output$ESS.L)

