#'
#'
#' Random-Walk Behaviour
#'
#'


# ==== Preamble ====

source("RWM Algorithm/Multivariate RWM Algorithm .R")

# ==== Stretched Gaussian RWM ====

Stretched.RWM = Multivariate.RWM(no.its = 1000, lambda = 1.6, 
                                 target = "Stretched.Gaussian", proposal = "Std.Gaussian",
                                 x0 = rmvn(1, mu = rep(0,2), sigma = diag(c(1,10), 2)), prop.V = diag(c(1,10), 2), burn.in = 0)

par(mfrow = c(2,2))

plot(Stretched.RWM$sample[,1], type = 'l')

plot(Stretched.RWM$sample[,2], type = 'l')

acf(Stretched.RWM$sample[,1])
acf(Stretched.RWM$sample[,2])

RWM = Multivariate.RWM(no.its = 1000, lambda = 1.4, 
                      target = "Std.Gaussian", proposal = "Std.Gaussian",
                      x0 = rmvn(1, mu = rep(0,2), sigma = diag(1, 2)), prop.V = diag(1, 2), burn.in = 0)

par(mfrow = c(2,2))

plot(RWM$sample[,1], type = 'l')

plot(RWM$sample[,2], type = 'l')

acf(RWM$sample[,1])
acf(RWM$sample[,2])