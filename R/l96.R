source("ode.R")


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

l96.test <- function(a=1e-5) {
  ns <- 40
  F <- 8
  dt <- 0.05
  nstep.spinup <- 1000
  seed <- 514
  set.seed(seed)
  
  x0 <- ode.fom(l96, rnorm(ns), nstep.spinup, dt, F)
  x <- ode.fom(l96, x0, 1, dt, F)
  
  dx0 <- a * rnorm(ns)
  x1 <- ode.fom(l96, x0 + dx0, 1, dt, F)
  dx <- ode.tlm(l96, l96.tl, x0, dx0, 1, dt, F)
  e <- sqrt(sum((x1 - x - dx)^2))
  cat("TLM:", "a=", a, "l2=", e, "\n")
  
  xa <- ode.adm(l96, l96.ad, x0, dx, 1, dt, F)
  tdxdx <- t(dx) %*% dx
  dx0xa <- t(dx0) %*% xa
  cat("ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa, "\n")
}