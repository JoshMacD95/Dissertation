#'
#' Further Topics in MCMC
#' Comparing Euler, Modified Euler and Leapfrog 
#'

# ===== Run Functions File ====

setwd('H:/MSc Statistics/Dissertation/Hamiltonian Monte Carlo')
#setwd('/Volumes/macdona8/MSc Statistics/Dissertation/Hamiltonian Monte Carlo')
source('Numerical Methods for systems of Diff Eqns.R')

t.seq = seq(0, 2*pi, by = 0.01)

q.seq = sapply(t.seq, FUN = q_t, q0 = 1, m = 1)

rho.seq = sapply(t.seq, FUN = rho_t, q0 = 1, m = 1)

# ==== Euler Method ====

Euler.test = Euler(q0 = 1, rho0 = 0, m = 1, timestep = 0.3, obs.time = 4*pi)

H.euler = c()

for(i in 1:nrow(Euler.test)){
  H.euler[i] = H(m = 1, q = Euler.test[i,1], rho = Euler.test[i,2])
}

# ==== Modified Euler Method ====

Euler.mod.test = Euler.mod(q0 = 1, rho0 = 0, m = 1, timestep = 0.3, obs.time = 4*pi)

H.euler.mod = c()

for(i in 1:nrow(Euler.mod.test)){
  H.euler.mod[i] = H(m = 1, q = Euler.mod.test[i,1], rho = Euler.mod.test[i,2])
}

# ==== Leapfrog Method ====

leapfrog.test = leapfrog(q0 = 1, rho0 = 0, m = 1, timestep = 0.3, obs.time = 4*pi)

H.leapfrog = c()

for(i in 1:nrow(leapfrog.test)){
  H.leapfrog[i] = H(m = 1, q = leapfrog.test[i,1], rho =leapfrog.test[i,2])
}

# ==== Comparing Results ====

par(mfrow = c(1,3))

plot(Euler.test[,1], Euler.test[,2], type = 'l')
lines(q.seq, rho.seq, type = 'l', col = 'blue')


plot(Euler.mod.test[,1], Euler.mod.test[,2], type = 'l')
lines(q.seq, rho.seq, type = 'l', col = 'blue')

plot(leapfrog.test[,1], leapfrog.test[,2], type = 'l')
lines(q.seq, rho.seq, type = 'l', col = 'blue')

par(mfrow = c(1,1))
plot(x = seq(1, length(H.euler)), y = rep(H.euler[1], length(H.euler)), type = 'l', col = 'green', lty = 2, ylim = c(-max(H.euler), max(H.euler)))
lines(H.euler, type = 'l', col = 'red')
lines(H.leapfrog, type = 'l', col = 'blue')
lines(H.euler.mod, type = 'l', col = 'purple')
