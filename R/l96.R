source("rk4.R")


l96 <- function(t, x, F=8) {
  n <- length(x)
  (x[c(2:n, 1)] - x[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] - x + F
}

l96.tl <- function(t, x, dx, ...){
  n <- length(x)
  (dx[c(2:n, 1)] - dx[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] +
    (x[c(2:n, 1)]- x[c(n-1, n, 1:(n-2))]) * dx[c(n, 1:(n-1))] - dx
}

l96.ad <- function(t, x, dxa, ...){
  n <- length(dxa)
  -x[c(2:n, 1)] * dxa[c(3:n, 1, 2)] + 
      (x[c(3:n, 1, 2)] - x[c(n, 1:(n-1))]) * dxa[c(2:n, 1)] +
      x[c(n-1, n, 1:(n-2))] * dxa[c(n, 1:(n-1))] - dxa
}

l96.fom <- function(x, nstep, dt, F) {
  t <- 0
  for (i in 1:nstep) {
    x <- rk4(l96, t, x, dt, F)
    t <- t + dt
  }
  x
}

l96.tlm <- function(x, dx, nstep, dt, F) {
  t <- 0
  for (i in 1:nstep) {
    dx <- rk4.tl(l96, l96.tl, t, x, dx, dt, F)
    t <- t + dt
  }
  dx
}

l96.adm <- function(x, xa, nstep, dt, F) {
  t <- nstep
  for (i in nstep:1){
    xa <- rk4.ad(l96, l96.ad, t, x, xa, dt, F)
    t <- t - dt
  }
  xa
}

test.l96 <- function(a=1e-5) {
  ns <- 40
  F <- 8
  dt <- 0.05
  nstep.spinup <- 1000
  seed <- 514
  set.seed(seed)
  
  x0 <- l96.fom(rnorm(ns), nstep.spinup, dt, F)
#  x0 <- rep(1, ns)
  x <- l96.fom(x0, 1, dt, F)
  
  dx0 <- a * rnorm(ns)
#  dx0 <- rep(a, ns)
  x1 <- l96.fom(x0 + dx0, 1, dt, F)
  dx <- l96.tlm(x0, dx0, 1, dt, F)
  e <- sqrt(sum((x1 - x - dx)^2))
  cat("TLM:", "a=", a, "l2=", e, "\n")
  
  xa <- l96.adm(x0, dx, 1, dt, F)
  tdxdx <- t(dx) %*% dx
  dx0xa <- t(dx0) %*% xa
  cat("ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa, "\n")
}