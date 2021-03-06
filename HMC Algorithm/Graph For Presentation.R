#'
#' Graphs for Presentation
#'
#'

# ==== Preamble ====

source("Numerical Methods/Numerical Methods for systems of Diff Eqns.R")
source("HMC Algorithm/Multivariate HMC Algorithm.R")

# ==== Creating Data ====
# Creating Density Values
x = seq(-4, 4, length = 10000)
pi_x = dnorm(x)
log_pi_x = log(pi_x)
neg_log_pi_x = -log(pi_x)
density.data = data.frame(x, pi_x, log_pi_x, neg_log_pi_x)


# Get positions for Hamiltonian Dynamics
leapfrog.pts = numerical.method(x0 = rnorm(1), rho0 = 3, m = 1, L = 100, obs.time = 2*pi, target ="1D.Gaussian",
                                K = squared.kinetic, method = leapfrog, final = FALSE)$x.values

# ==== Plots ====

# Exhibiting Hamiltonian Dynamics
for(i in 1:length(leapfrog.pts)){
  suffix = c(rep(0, nchar(length(leapfrog.pts)) - nchar(i)), i, "", "", "")
  plotname = paste('MCMC using Hamiltonian Dynamics Presentation/Animations/First_Normal_Potential', suffix[1], suffix[2], suffix[3], suffix[4], ".png" , sep = "" )
  png(plotname, width = 680, height = 480, units = 'px')
  par(mfrow = c(1,1))
  #plot(x, pi_x, type = 'l', col = 'blue', ylab = expression(pi(x)), main = "N(0,1) Target")
  #plot(x, log_pi_x, type = 'l', col = 'blue', xlab = expression(x), ylab = expression(log(pi(x))), main = "log(N(0,1))")
  plot(x, neg_log_pi_x, type = 'l', ylab = expression(U(x)), col = 'blue', main = " Potential Surface U(x)")
  points(leapfrog.pts[i], U(leapfrog.pts[i], target = "1D.Gaussian"), col = 'tomato2', cex = 2, pch =16)
  dev.off()
}

# Plotting RWM Optimal Acceptance Rate graphs

RWM.output2 = read.csv("Output Data/Optimal_RWMDimension2.csv")
RWM.output4 = read.csv("Output Data/Optimal_RWMDimension4.csv")
RWM.output8 = read.csv("Output Data/Optimal_RWMDimension8.csv")
RWM.output12 = read.csv("Output Data/Optimal_RWMDimension12_2.csv")

plot(RWM.output2$AR, RWM.output2$ESS.sec, main = "RWM Algorithm, Multivariate Normal Target Dimension 2",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate')

abline(v = 0.234, col = 'red', lty = 2)

plot(RWM.output4$AR, RWM.output4$ESS.sec, main = "RWM Algorithm, Multivariate Normal Target Dimension 4",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate')

abline(v = 0.234, col = 'red', lty = 2)


plot(RWM.output8$AR, RWM.output8$ESS.sec, main = "RWM Algorithm, Multivariate Normal Target Dimension 8",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate', ylim = c(0,700))

abline(v = 0.234, col = 'red', lty = 2)

plot(RWM.output12$AR, RWM.output12$ESS.sec, main = "RWM Algorithm, Multivariate Normal Target Dimension 12",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate', ylim = c(0,300))

abline(v = 0.234, col = 'red', lty = 2)

HMC.output1 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (2, 4, 10, 50, 100).csv")
HMC.output2 = HMC.output1[HMC.output1$Dimension==2,]
HMC.output100 = HMC.output1[HMC.output1$Dimension==100,]
HMC.output500 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (500)2.csv")
HMC.output1000 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension(1000).csv")

# ==== ESS vs. L ====

mean.ESS = c()
mean.ESS.L = c()
for(i in 1:20){
  mean.ESS.L[i] = mean(HMC.output2$ESS.L[HMC.output100$L==i])
  mean.ESS[i] = mean(HMC.output2$ESS[HMC.output100$L==i])
}
par(mfrow = c(1,2))
plot(1:20, mean.ESS, ylab = "Mean ESS", xlab = "L", pch = 16, col ='light blue')
plot(1:20, mean.ESS.L, ylab = "Mean ESS/L", xlab = "L", ylim = c(0, 1300), pch = 16, col ='tomato2')
# ==== Optimal Acceptance Rate Graphs ====
plot(HMC.output100$AR, HMC.output100$ESS.L, main = "HMC Algorithm, Stetched Logistic Target Dimension 100",
     ylab = 'ESS/L', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)

plot(HMC.output500$AR, HMC.output500$ESS.L, main = "HMC Algorithm, Stretched Logistic Target Dimension 500",
     ylab = 'ESS/L', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)


plot(HMC.output1000$AR, HMC.output1000$ESS.L, main = "HMC Algorithm, Stretched Logistic Target Dimension 1000",
     ylab = 'ESS/L', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)



# Function which creates ggplot of Potential Surface and the current position
Univariate.potential.plot = function(dataset, position, target = "1D.Gaussian", target.colour = "blue", position.colour = "tomato2"){
  p = ggplot(data = dataset, aes(x = x, y = neg_log_pi_x)) + geom_line(color = target.colour) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab(expression(x)) + ylab(expression(U(x))) +
    geom_point(aes(position, U(position, target)), color = position.colour, size = 4) 
  print(p)
}

# Test Plots
#Univariate.potential.plot(density.data, leapfrog.pts[1], target.colour = "dark blue")

# ==== Plotting ====

# Create Blank Canvas for plots in animation
img = image_graph(600, 340, res = 96)

# Create plots which will be the frames of the animation
plots = lapply(leapfrog.pts, Univariate.potential.plot, dataset = density.data)

# Shuts down current device?
dev.off()

# ==== Creating Animation ====

# Create Animation from collection of plots in image
animation = image_animate(img, fps = 100)

# Preview the Animation
print(animation)

# Save Animation as a gif
image_write(animation, "Graphs/Animations/Ball Rolling on Potential Surface.gif")
