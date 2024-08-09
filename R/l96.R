source("rk4.R")


l96 <- function(t, x, F=8) {
  n <- length(x)
  (x[c(2:n, 1)] - x[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] - x + F
}

l96.tl <- function(t, x, dx){
  n <- length(x)
  (dx[c(2:n, 1)] - dx[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] +
    (x[c(2:n, 1)]- x[c(n-1, n, 1:(n-2))]) * dx[c(n, 1:(n-1))] - dx
}

l96.ad <- function(t, x, xa){
  n <- length(xa)
  xa <- xa - x[c(2:n, 1)] * xa[c(3:n, 1, 2)] + 
    (x[c(3:n, 1, 2)] - x[c(n, 1:(n-1))]) * xa +
    x[c(n-1, n, 1:(n-2))] * xa[c(n, 1:(n-1))] - xa
  xa = 0.0
  xa
}

l96.fom <- function(x, nstep, dt, F) {
  t <- 0
  for (i in 1:nstep) {
    x <- rk4(l96, t, x, dt, F)
    t <- t + dt
  }
  x
}

l96.tlm <- function(x, dx, nstep, dt) {
  t <- 0
  for (i in 1:nstep) {
    dx <- rk4.tl(l96, l96.tl, t, x, dx, dt)
    t <- t + dt
  }
  dx
}

l96.adm <- function(x, xa, nstep, dt) {
  t <- nstep
  for (i in nstep:1){
    xa <- rk4.ad(l96, l96.ad, t, x, xa, dt)
    t <- t - dt
  }
  xa
}

test.l96 <- function() {
  ns <- 40
  F <- 8
  dt <- 0.05
  nstep.spinup <- 1000
  a <- 1e-5
  seed <- 514
  set.seed(seed)
  
  x0 <- l96.fom(rnorm(ns), nstep.spinup, dt, F)
  x <- l96.fom(x0, 1, dt, F)
  
  dx <- a * rnorm(ns)
  x1 <- l96.fom(x0 + dx, 1, dt, F)
  dx <- l96.tlm(x0, dx, 1, dt)
  e <- sqrt(sum((x1 - x - dx)^2))
  cat("TLM check:", "a=", a, "error=", e, "\n")
  
  xa <- l96.adm(x0, dx, 1, dt)
  sprintf("ADJ check: %e", (t(dx) %*% dx)-(t(dx) %*% xa))
}