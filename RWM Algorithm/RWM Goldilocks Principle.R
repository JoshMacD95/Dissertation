#'
#' Showing the Goldilocks Principle
#' 
#'

# ==== Preamble ==== 

source("RWM Algorithm/Multivariate RWM Algorithm .R")


# ==== Simulation ====

# Small Proposal Parameter

lambda.small = 0.1

small.RWM = Multivariate.RWM(no.its = 10000, lambda = lambda.small, 
                             target = "Std.Gaussian", proposal = "Std.Gaussian",
                             x0 = rnorm(1), burn.in = 0)

plot(small.RWM$sample, type = 'l')

# Large Proposal Parameter 

lambda.large = 100


large.RWM = Multivariate.RWM(no.its = 10000, lambda = lambda.large, 
                             target = "Std.Gaussian", proposal = "Std.Gaussian",
                             x0 = rnorm(1), burn.in = 0)

plot(large.RWM$sample, type = 'l')



# Medium Proposal Parameter

lambda.medium = 2.8


medium.RWM = Multivariate.RWM(no.its = 10000, lambda = lambda.medium, 
                             target = "Std.Gaussian", proposal = "Std.Gaussian",
                             x0 = rnorm(1), burn.in = 0)



# ==== Graphs ====

par(mfrow= c(2,3))

# Trace Plots
plot(small.RWM$sample, xlab = expression(i), ylab = expression(x[i]), type = 'l')
plot(large.RWM$sample, xlab = expression(i), ylab = expression(x[i]), type = 'l')
plot(medium.RWM$sample, xlab = expression(i), ylab = expression(x[i]), type = 'l')

# ACF Plots
acf(small.RWM$sample, main = '')
acf(large.RWM$sample, main = '')
acf(medium.RWM$sample, main = '')


plot(medium.RWM$sample, type = 'l')

qqplot(rnorm(length(medium.RWM$sample)),medium.RWM$sample)

hist(medium.RWM$sample, freq = F)
lines(seq(-4, 4, length = 10000), dnorm(seq(-4, 4, length = 10000)), col = 'blue')





