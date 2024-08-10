source("l96.R")

seed <- 514
set.seed(seed)

ns <- 40
F <- 8

nw <- 4
nt <- 1000
nc <- 250
nt <- nc * nw
ni <- 100

r <- 4
b <- 1
dt <- 0.05
a <- 0.1
g.break <- 1e-8

calc.cost <- function(dx, b, dy, r) {
  cost <- 0.5 * t(dx) %*% dx / b
  for (i in 1:ncol(dy)) {
    cost <- cost + 0.5 * t(dy[, i]) %*% dy[, i] / r
  }
  as.numeric(cost)
}

xt <- l96.fom(rnorm(ns), nt, dt, F)
x0 <- l96.fom(rnorm(ns), nt, dt, F)

yo <- matrix(rep(0, ns*nw), ncol=nw)
xb <- matrix(rep(0, ns*nw), ncol=nw)
d <- matrix(rep(0, ns*nw), ncol=nw)
l2 <- rep(0, nc)

for (k in 1:nc){
  cat("cycle=", k, "\n")
  xt0 = xt
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
      ad <- l96.adm(xb[, j], ad, 1, dt, F) + d[, j] / r
    }
    xb[, 1] <- xb[, 1] - a * ad
    cost <- calc.cost(xb[, 1] - x0, b, d, r)
    if (i == 1) {
      cost.old <- cost
    }
    l2[k] <- sqrt(mean((xb[, 1] - xt0)^2))
    gnorm <- t(ad) %*% ad
    cat(sprintf("%d J=%e Jold=%e g=%e l2=%e\n", i, cost, cost.old, gnorm, l2[k]))
    if (gnorm < g.break) {
      break
    }
    if (cost > cost.old) {
      xb[, 1] <- xb[, 1] + a * ad
      e[k] <- sqrt(mean((xb[, 1] - xt0)^2))
      break
    }
    cost.old <- cost
  }
  x0 <- l96.fom(xb[, 1], nw + 1, dt, F)
}

title <- paste("4DVar", "step=", a)
plot(seq(1, nt, by=nw), l2, type="l", main=title,
     xlab="time", ylab="l2", ylim=c(0, 3))
abline(h=1, lt=2)
abline(h=mean(l2), lt=3)
legend("topleft", legend=c("analysis", paste("mean", format(mean(e), digits=3)), "observation"),
       lty=c(1, 3, 2))
