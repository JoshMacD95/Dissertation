#''
#'
#' Optimal Acceptance Rate Graphs
#'
#'

HMC.output.lower = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (2, 4, 10, 50, 100).csv")

HMC.output.100 = HMC.output.lower[HMC.output.lower$Dimension == 100,]
HMC.output.500 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension (500).csv")
HMC.output.1000 = read.csv("Output Data/Logistic_HMC_With_Increasing_Dimension(1000).csv")



plot(HMC.output.100$AR, HMC.output.100$ESS.L, main = "HMC Algorithm, Stretched Logistic Target Dimension 100",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)

plot(HMC.output.500$AR, HMC.output.500$ESS.L, main = "HMC Algorithm, Stretched Logistic Target Dimension 500",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)

plot(HMC.output.1000$AR, HMC.output.1000$ESS.L, main = "HMC Algorithm, Stretched Logistic Target Dimension 1000",
     ylab = 'ESS/sec', xlab = 'Acceptance Rate')

abline(v = 0.65, col = 'red', lty = 2)



median.ESS.L = c()
UQ.ESS.L = c()
LQ.ESS.L = c()
median.AR = c()
for(i in 1:20){
  median.ESS.L[i] = median(output.1000$ESS.L[output.1000$L == i])
  UQ.ESS.L[i] = quantile(output.1000$ESS.L[output.1000$L == i], probs = 0.75) 
  LQ.ESS.L[i] = quantile(output.1000$ESS.L[output.1000$L == i], probs = 0.25) 
  median.AR[i] = median(output.1000$AR[output.1000$L == i])
}

plot(median.AR, median.ESS.L, ylim = c(0, max(UQ.ESS.L)), xlim = c(0,1))

points(median.AR, LQ.ESS.L, col = 'red')
points(median.AR, UQ.ESS.L, col = 'red')

median.ESS.L.100 = c()
UQ.ESS.L.100 = c()
LQ.ESS.L.100 = c()
median.AR.100 = c()
for(i in 1:20){
  median.ESS.L.100[i] = median(output.100$ESS.L[output.100$L == i])
  UQ.ESS.L.100[i] = quantile(output.100$ESS.L[output.100$L == i], probs = 0.75) 
  LQ.ESS.L.100[i] = quantile(output.100$ESS.L[output.100$L == i], probs = 0.25) 
  median.AR.100[i] = median(output.100$AR[output.100$L == i])
}


plot(median.AR.100, median.ESS.L.100, ylim = c(0, max(UQ.ESS.L.100)), xlim = c(0,1))

points(median.AR.100, LQ.ESS.L.100, col = 'red')
points(median.AR.100, UQ.ESS.L.100, col = 'red')

abline(v = 0.8, col = 'blue', lty = 2)
L.list = 1:20
d.list = c(2, 4, 10, 50, 100, 500, 1000)

d = c()
L = c()
for(i in  1:length(L.list)){
  for(j in 1:length(d.list)){
    
    mean.ESS.L = mean
    
  }
}


write.csv(output, file = "Optimal Acceptance Rate Results")

d = c(2, 4, 10, 50, 100, 500, 1000)
par(mfrow=c(1,1))
for(i in 1:length(d)){
  output.d = output[output$Dimension == d[i],]
  plot(output.d$AR, output.d$ESS.L, 
       ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = paste("Efficiency of Logisitic Sampling in Dimension", d[i], sep = " "),
       col = 1, 
       ylim = c(0, max(output.d$ESS.L)), xlim = c(0,1))
  abline(v = output.d$AR[which.max(output.d$ESS.L)], col = 'red', lty = 2)
}
plot(output$AR[output$Dimension==d[i]], output$ESS.L[output$Dimension==d[i]], 
     ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = paste("Efficiency of Logisitic Sampling in Dimension", d[i], sep = ""),
     col = 1, 
     ylim = c(0, max(output$ESS.L[output$Dimension==2])), xlim = c(0,1))
plot(output$AR[output$Dimension==4], output$ESS.L[output$Dimension==4], 
     ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = "Efficiency of Logisitic Sampling in Dimension 4",
     col = 1, 
     ylim = c(0, max(output$ESS.L[output$Dimension==4])), xlim = c(0,1))
plot(output$AR[output$Dimension==10], output$ESS.L[output$Dimension==10], 
     ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = "Efficiency of Logisitic Sampling in Dimension 10",
     col = 1, 
     ylim = c(0, max(output$ESS.L[output$Dimension==10])), xlim = c(0,1))
plot(output$AR[output$Dimension==50], output$ESS.L[output$Dimension==50], 
     ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = "Efficiency of Logisitic Sampling in Dimension 50",
     col = 1, 
     ylim = c(0, max(output$ESS.L[output$Dimension==50])), xlim = c(0,1))
plot(output$AR[output$Dimension==100], output$ESS.L[output$Dimension==100], 
     ylab = "ESS per Leapfrog Step", xlab = "Acceptance Rate", main = "Efficiency of Logisitic Sampling in Dimension 100",
     col = 1, 
     ylim = c(0, max(output$ESS.L[output$Dimension==100])), xlim = c(0,1))



