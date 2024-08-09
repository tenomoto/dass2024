source("l96.R")

seed <- 514
set.seed(seed)

ns <- 40
F <- 8

nw <- 4
nt.spinup <- 1000
nt <- 1000
nc <- nt %/% nw
ni <- 100

r <- 4
b <- 1
dt <- 0.05
a <- 0.8
g.break <- 1e-8

calc.cost <- function(dx, b, dy, r) {
  cost <- 0.5 * t(dx) %*% dx / b
  for (i in 1:ncol(d)) {
    cost <- cost + 0.5 * t(dy[, j]) %*% dy[, j] / r
  }
  as.numeric(cost)
}

xt <- l96.fom(rnorm(ns), nt.spinup, dt, F)
x0 <- l96.fom(rnorm(ns), nt.spinup, dt, F)

yo <- matrix(rep(0, ns*nw), ncol=nw)
xb <- matrix(rep(0, ns*nw), ncol=nw)
d <- matrix(rep(0, ns*nw), ncol=nw)
e <- rep(0, nc)

for (k in 1:nc){
  cat("cycle=", k, "\n")
  for (j in 1:nw){
    yo[, j] <- rnorm(ns, xt, sqrt(r))
    xt <- l96.fom(xt, 1, dt, F)
  }
  xb[, 1] <- x0
  for (i in 1:ni) {
    for (j in 2:nw){
      xb[, j] <- l96.fom(xb[, j-1], 1, dt, F)
    }
    d = xb - yo
    ad <- rep(0, ns)
    for (j in nw:1){
      ad <- l96.adm(xb[, j], ad, 1, dt) + d[, j] / r
    }
    xb[, 1] <- xb[, 1] - a * ad
    cost <- calc.cost(xb[, 1] - x0, b, d, r)
    if (i == 1) {
      cost.old <- cost
    }
    e[k] <- sqrt(mean(xb[, 1] - xt)^2)
    gnorm <- t(ad) %*% ad
    cat(sprintf("%d J=%e g=%e l2=%e\n", i, cost, gnorm, e[k]))
    if (gnorm < g.break) {
      break
    }
    if (cost > cost.old) {
      xb[, 1] <- xb[, 1] + a * ad
      break
    }
    cost <- cost.old
  }
  x0 <- l96.fom(xb[, 1], nw + 1, dt, F)
}

title <- paste("4DVar", "step=", a)
plot(seq(1, nt, by=nw), e, type="l", main=title,
     xlab="time", ylab="error", ylim=c(0, 3))
abline(h=1, lt=2)
abline(h=mean(e), lt=3)
legend("topleft", legend=c("analysis", paste("mean", format(mean(e), digits=3)), "observation"),
       lty=c(1, 3, 2))