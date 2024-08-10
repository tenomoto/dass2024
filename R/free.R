source("l96.R")
source("ode.R")

nj <- 40
nstep <- 1001
x.hist <- matrix(0, nj, nstep)
x <- rnorm(nj)
F <- 3.85
dt <- 0.05
t <- 0.0
for (i in 1:nstep-1) {
  x <- rk4(l96, t, x, dt, F)
  x.hist[,i] <- x
  t <- t + dt
}
t <- seq(0, nstep*dt, length.out=nstep)
filled.contour(1:nj, t, x.hist, nlevel=11,
               ylim=rev(range(t)), xlab="j", ylab="time")