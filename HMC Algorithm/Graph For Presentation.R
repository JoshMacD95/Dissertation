#'
#' Graphs for Presentation
#'
#'

# ==== Preamble ====

source("Numerical Methods/Numerical Methods for systems of Diff Eqns.R")

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

for(i in 1:length(leapfrog.pts)){
  suffix = c(rep(0, nchar(length(leapfrog.pts)) - nchar(i)), i, "", "", "")
  plotname = paste('MCMC using Hamiltonian Dynamics Presentation/Animations/Normal_Potential', suffix[1], suffix[2], suffix[3], suffix[4], ".png" , sep = "" )
  png(plotname, width = 680, height = 480, units = 'px')
  par(mfrow = c(1,3))
  plot(x, pi_x, type = 'l', col = 'blue', ylab = expression(pi(x)), main = "N(0,1) Target")
  plot(x, log_pi_x, type = 'l', col = 'blue', xlab = expression(x), ylab = expression(log(pi(x))), main = "log(N(0,1))")
  plot(x, neg_log_pi_x, type = 'l', ylab = expression(U(x)), col = 'blue', main = " Potential Surface U(x)")
  points(leapfrog.pts[i], U(leapfrog.pts[i], target = "1D.Gaussian"), col = 'tomato2', cex = 2, pch =16)
  dev.off()
}























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
