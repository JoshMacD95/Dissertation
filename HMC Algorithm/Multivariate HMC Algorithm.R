#'
#' Further Topics in MCMC
#' Multivariate HMC Algortihm
#'

# ==== Preamble ====
library(coda)
library(mvnfast)
source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')

# ==== Multivariate Hamiltonian Algorithm ====

Multivariate.HMC = function(q0, m, L, obs.time, 
                            U, grad.U, K, method = leapfrog, 
                            no.its, burn.in, output.type = "all"){
  
  time.start = Sys.time()
  
  d = length(q0)
  
  if(length(m) == 1){
    m = rep(m, d)
  }
  if(length(m) != d){
    return(print("Mass vector not same dimension as Position vector"))
  }
  
  # Initial Potential Energy
  U0 = U(q0)
  
  q.curr = q0
  U.curr = U0
  
  no.accept = 0
  draws = matrix(data = 0, nrow = no.its + 1, ncol = d + 1)
  
  draws[1,] = c(q.curr, U.curr)
  
  
  for(i in 2:(no.its+1)){
    
    # Momentum Step
    rho.curr = rmvn(1, mu = rep(0,d), sigma = diag(m))
    K.curr = K(rho.curr, m)
    
    # Hamiltonian Dynamics 
    ham.dyn = numerical.method(q.curr, rho.curr, m, L, obs.time, grad.U, method, final = TRUE)
    q.prop = ham.dyn[,1]
    rho.prop = - ham.dyn[,2]
    
    U.prop = U(q.prop)
    K.prop = K(rho.prop, m)
    
    # Accept-recject Step
    
    if(runif(1) < exp(- U.prop - K.prop + U.curr + K.curr)){
      
      # Update position and -log(posteriors)
      q.curr = q.prop
      U.curr = U.prop
      rho.curr = rho.prop
      
      no.accept = no.accept + 1
    }
    
    draws[i,] = c(q.curr, U.curr)
  }
  draws = draws[-(1:burn.in),]
  time.finish = Sys.time()
  time.taken = time.finish - time.start
  # Calculating Effective Sample Size
  ESS = effectiveSize(draws[,1:d])
  ESS.per.sec = ESS/as.numeric(time.taken)
  # Calculating Acceptance Rate
  accept.rate = no.accept/no.its
  
  # Important Output
  output(U, K, ESS, accept.rate, q0 = q0, no.its)
  if(output.type == "all"){
    return(list(sample = draws[,1:d], log.target = -draws[,d+1], ESS = ESS, ESS.per.sec =ESS.per.sec))
  }
  if(output.type == "ESS"){
    return(ESS)
  }
  if(output.type == "ESS/sec"){
    return(ESS.per.sec)
  }
}

Test.MHMC = Multivariate.HMC(q0 = c(0,0), L = 100, obs.time = pi/2, m = 1, 
                 U = MultGauss, grad.U = grad.MultGauss, K = squared.kinetic,
                 no.its = 10000, burn.in = 1000, output.type = "all")

Test.MHMC
par(mfrow= c(1,2))
plot(Test.MHMC$sample[,1], type = 'l')

plot(Test.MHMC$sample[,2], type = 'l')

x_1 = seq(from = -3, to = 3, length = 1000)
x_2 = seq(from = -3, to = 3, length = 1000)
mvn.values = matrix(data = NA, nrow = 1000, ncol = 1000)

for(i in 1:1000){
  for(j in 1:1000){
    mvn.values[i,j] = dmvn(X = c(x_1[i], x_2[j]), mu = rep(0, 2), sigma = diag(rep(1, 2)))
  }
}

contour(x_1, x_2, mvn.values)
points(Test.MHMC$sample[,1], Test.MHMC$sample[,2])
m = 1

d = length(q0)

if(length(m) == 1){
  m = rep(m, d)
}
if(length(m) != d){
  return(print("Mass vector not same dimension as Position vector"))
}

rho.curr = rmvn(1, mu = rep(0,d), sigma = diag(m))
# Initial Potential Energy
U0 = MultGauss(q=c(0,0))

q.curr = q0
U.curr = U0

no.accept = 0
draws = matrix(data = 0, nrow = no.its + 1, ncol = d + 1)

draws[1,] = c(q.curr, U.curr)


