source("l96.R")
source("step.R")
source("config_l96.R")

nc <- nt %/% nw
calc.cost <- function(dx, b, dy, r) {
  cost <- 0.5 * t(dx) %*% dx / b
  for (j in 1:ncol(dy)) {
    cost <- cost + 0.5 * t(dy[, j]) %*% dy[, j] / r
  }
  as.numeric(cost)
}

con.true <- file(xt_fname, "rb")
con.obs <- file(yo_fname, "rb")
con <- file(xf_fname, "rb")
x0 <- readBin(con, "numeric", ns)
close(con)

xb <- matrix(rep(0, ns*nw), ncol=nw)
l2 <- rep(0, nc)

for (k in 1:nc){
  cat("cycle=", k, "\n")
  xt0 <- readBin(con.true, "numeric", ns)
  for (j in 2:nw){
      xt <- readBin(con.true, "numeric", ns)
  }
  yo <- matrix(readBin(con.obs, "numeric", ns*nw), nrow=ns)
  xb[, 1] <- x0
  for (i in 1:ni) {
    for (j in 2:nw){
      xb[, j] <- step.fom(l96, xb[, j-1], 1, dt, F)
    }
    d = xb - yo
    ad <- rep(0, ns)
    for (j in nw:1){
      ad <- step.adm(l96, l96.ad, xb[, j], ad, 1, dt, F) + d[, j] / r
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
      l2[k] <- sqrt(mean((xb[, 1] - xt0)^2))
      break
    }
    cost.old <- cost
  }
  x0 <- step.fom(l96, xb[, 1], nw + 1, dt, F)
}
close(con.true)
close(con.obs)

con.out <- file("l96_var.dat", "wb")
writeBin(l2, con.out)
close(con.out)
