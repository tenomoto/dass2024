source("rk4.R")

step.fom <- function(f, x, nstep, dt, ...) {
  t <- 0
  for (i in 1:nstep) {
    x <- rk4(f, t, x, dt, ...)
    t <- t + dt
  }
  x
}

step.tlm <- function(f, df, x, dx, nstep, dt, ...) {
  t <- 0
  for (i in 1:nstep) {
    dx <- rk4.tl(f, df, t, x, dx, dt, ...)
    t <- t + dt
  }
  dx
}

step.adm <- function(f, fa, x, xa, nstep, dt, ...) {
  t <- nstep
  for (i in nstep:1){
    xa <- rk4.ad(f, fa, t, x, xa, dt, ...)
    t <- t - dt
  }
  xa
}
