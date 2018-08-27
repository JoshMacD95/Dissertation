#'
#'
#'
#' Optimal Acceptance Rate Graphs (using ggplot2)
#'
#'

# ==== Preamble ====
library(ggplot2)
library(RColorBrewer)
library(grid)
# ==== Datasets ====

# == RWM Optimal Acceptance Rate Output ==

RWM.output1 = read.csv("Output Data/Optimal_RWMDimension1.csv")
RWM.output2 = read.csv("Output Data/Optimal_RWMDimension2.csv")
RWM.output4 = read.csv("Output Data/Optimal_RWMDimension4.csv")
RWM.output8 = read.csv("Output Data/Optimal_RWMDimension8.csv")
RWM.output12 = read.csv("Output Data/Optimal_RWMDimension12_2.csv")

# == HMC Optimal Acceptance Rate Output ==

HMC.output1 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (2, 4, 10, 50, 100).csv")
HMC.output2 = HMC.output1[HMC.output1$Dimension==2,]
HMC.output100 = HMC.output1[HMC.output1$Dimension==100,]
HMC.output500 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (500)2.csv")
HMC.output1000 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension(1000).csv")


# ==== Plots ====

# == 1D RWM Optimal ==


# == Specify Data and Variables
data = ggplot(data = RWM.output1, aes(x = AR, y = ESS.sec))

# == Specify type of plot and variables to be grouped or coloured by
points = geom_point(aes(col= lambda)) 

# == Specify x and y axis labels, limits and ticks etc.
x.axis = scale_x_continuous(name =expression(alpha[rate]), 
                          limits=c(0,1))
y.axis = scale_y_continuous(name = "ESS/sec")

# == Specify a legend/colourscale for colour coding + parameters
scale = scale_color_continuous(minor_breaks = NULL, low = "light blue", high = "dark blue", guide = "colourbar", name = expression(lambda), limits = c(0,10)) 
scale.parameters = guides(colour = guide_colorbar(title.position = "left", title.hjust = 1, ticks = FALSE))

# == Specify Paramters of the plot grid
grid.parameters = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = c(0.8, 0.8) ) 

# == Add line at optimal acceptance rate 
opt.acc.line = geom_vline(xintercept = 0.44, color = 'red', linetype = 'dashed')
# == Plot graph
data + points + x.axis + y.axis + legend + grid.parameters + scale + scale.parameters + opt.acc.line

# == 2, 4, 8, 12D RWM Optimal ==

# == Specify Data and Variables
data2 = ggplot(data = RWM.output2, aes(x = AR, y = ESS.sec))
data4 = ggplot(data = RWM.output4, aes(x = AR, y = ESS.sec))
data8 = ggplot(data = RWM.output8, aes(x = AR, y = ESS.sec))
data8 = ggplot(data = RWM.output12, aes(x = AR, y = ESS.sec))
# == Specify type of plot and variables to be grouped or coloured by
points = geom_point(aes(col = lambda)) 

# == Specify x and y axis labels, limits and ticks etc.
x.axis = scale_x_continuous(name =expression(alpha[rate]), 
                            limits=c(0,1))
y.axis = scale_y_continuous(name = "ESS/sec")

# == Specify a legend/colourscale for colour coding + parameters
scale12 = scale_color_continuous(minor_breaks = NULL, low = "light blue", high = "dark blue", guide = "colourbar", name = expression(lambda), limits = c(0,10)) 
scale.parameters12 = guides(colour = guide_colorbar(title.position = "left", title.hjust = 1, ticks = FALSE))

# == Specify Paramters of the plot grid
grid.parameters = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = c(0.8, 0.8) ) 

# == Add line at optimal acceptance rate 
opt.acc.line = geom_vline(xintercept = 0.234, color = 'red', linetype = 'dashed')

RWM2 = data2 + points + x.axis + y.axis + grid.parameters + opt.acc.line
RWM4 = data4 + points + x.axis + y.axis + grid.parameters + opt.acc.line
RWM4 = data8 + points + x.axis + y.axis + grid.parameters + opt.acc.line
RWM4 = data12 + points + x.axis + y.axis + grid.parameters + opt.acc.line
