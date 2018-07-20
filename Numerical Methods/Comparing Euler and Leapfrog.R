#'
#' Further Topics in MCMC
#' Comparing Euler, Modified Euler and Leapfrog 
#'

# ===== Run Functions File ====

source('Numerical Methods/Numerical Methods for systems of Diff Eqns.R')

t.seq = seq(0, 2*pi, by = 0.01)

x.seq = sapply(t.seq, FUN = x_t, x0 = 1, m = 1)

rho.seq = sapply(t.seq, FUN = rho_t, x0 = 1, m = 1)

# ==== Euler Method ====

Euler.test = numerical.method(x0 = 0, rho0 = 1, m = 1, L = 20, obs.time = 2*pi, target = "1D.Gaussian", method = Euler, final = FALSE)

H.euler = c()

for(i in 1:nrow(Euler.test)){
  H.euler[i] = U(x = Euler.test[1,1], target = "1D.Gaussian") + squared.kinetic(rho = Euler.test[i,2], m = 1)
}

# ==== Modified Euler Method ====

Euler.mod.test = numerical.method(x0 = 0, rho0 = 1, m = 1, L = 20, obs.time = 2*pi, target = "1D.Gaussian", method = Euler.mod, final = FALSE)

H.euler.mod = c()

for(i in 1:nrow(Euler.mod.test)){
  H.euler.mod[i] = U(x = Euler.mod.test[i,1], target = "1D.Gaussian") + squared.kinetic(rho = Euler.mod.test[i,2], m = 1)
}

# ==== Leapfrog Method ====

leapfrog.test = numerical.method(x0 = 0, rho0 = 1, m = 1, L = 20, obs.time = 2*pi, target = "1D.Gaussian", method = leapfrog, final = FALSE)
H.leapfrog = c()

for(i in 1:nrow(leapfrog.test)){
  H.leapfrog[i] = U(x = leapfrog.test[i,1], target = "1D.Gaussian") + squared.kinetic(rho = leapfrog.test[i,2], m = 1)
}

# ==== Comparing Results ====

par(mfrow = c(2,3))

plot(Euler.test[,1], Euler.test[,2], type = 'l',  ylab = 'Momentum', xlab = 'Position')
lines(x.seq, rho.seq, type = 'l', col = 'blue')


plot(Euler.mod.test[,1], Euler.mod.test[,2], type = 'l',  ylab = 'Momentum', xlab = 'Position')
lines(x.seq, rho.seq, type = 'l', col = 'blue')

plot(leapfrog.test[,1], leapfrog.test[,2], type = 'l',  ylab = 'Momentum', xlab = 'Position')
lines(x.seq, rho.seq, type = 'l', col = 'blue')


plot(x = seq(1, length(H.euler)), y = rep(H.euler[1], length(H.euler)), 
     type = 'l', col = 'green', lty = 2, 
     ylim= c(H.euler[1] - max(H.euler),H.euler[1] + max(H.euler)),
     xlab = "Step Number", ylab = "Hamiltonian"
     )
lines(H.euler, type = 'l', col = 'red')

plot(x = seq(1, length(H.euler.mod)), y = rep(H.euler.mod[1], length(H.euler.mod)), 
     type = 'l', col = 'green', lty = 2,
     xlab = "Step Number", ylab = "Hamiltonian"
     )
lines(H.euler.mod, type = 'l', col = 'purple')

plot(x = seq(1, length(H.leapfrog)), y = rep(H.leapfrog[1], length(H.leapfrog)), 
     type = 'l', col = 'green', lty = 2, 
     xlab = "Step Number", ylab = "Hamiltonian")
lines(H.leapfrog, type = 'l', col = 'blue')


