source("l63.R")
source("step.R")
source("config_l63.R")

seed <- 514
set.seed(seed)

calc.cost <- function(dx, b, dy, r) {
  cost <- 0.5 * t(dx) %*% dx / b
  for (j in 1:ncol(dy)) {
    cost <- cost + 0.5 * t(dy[, j]) %*% (dy[, j] / r)
  }
  as.numeric(cost)
}
options(digits = 5)
xt <- matrix(rep(0, ns * nt), nrow=ns)
xt[, 1] <- xt0
yo <- matrix(rep(0, ns * nobs), nrow=ns)
m <- 1
for (k in 2:nt) {
  xt[, k] <- step.fom(l63, xt[, k-1], 1, dt, p, r, b)
  if (k %% obs.int == 0) {
    yo[, m] <- rnorm(ns, xt[, k], sqrt(obs.r))
    m <- m + 1
  }
}

l2 <- rep(0, ni)
cost <- rep(0, ni)
xb <- matrix(rep(0, ns * nt), nrow=ns)
xb[, 1] <- xb0

for (i in 1:ni) {
  for (j in 2:nt) {
    xb[, j] <- step.fom(l63, xb[, j-1], 1, dt, p, r, b)    
  }
  d <- xb[, (1:nobs) * obs.int] - yo
  ad <- rep(0, ns)
  m = nobs
  for (j in nt:1){
    ad <- step.adm(l63, l63.ad, xb[, j], ad, 1, dt, p, r, b)
    if (j %% obs.int == 0) {
      ad <- ad + d[, m] / obs.r
      m <- m - 1
    }
  }
  xb[, 1] <- xb[, 1] - a * ad
  l2[i] <- sqrt(mean((xb[, 1] - xt0)^2))
  cost[i] <- calc.cost(xb[,1] - xb0, bgd.b, d, obs.r)
  cat(i, l2[i], cost[i], "\n")
}

tab10.new <- c('#5778a4', '#e49444', '#d1615d', '#85b6b2', '#6a9f58',
               '#e7ca60', '#a87c9f', '#f1a2a9', '#967662', '#b8b0ac')

title <- paste("L63 Var a=", a, " ni=", ni, " L2=", format(l2[ni], digits=3))
off <- c(0, 40, 60)
t <- (1:nt) * dt
to <- ((1:nobs) * obs.int) * dt
xt <- xt + off
xb <- xb + off
yo <- yo + off
plot(1, type = "n", xlab="time", ylab="state", main=title,
     xlim=c(min(t), max(t)), ylim=c(-20, 150))
for (i in 1:ns) {
  lines(t, xt[i,], lwd=2, lty=2, col=tab10.new[i])
  lines(t, xb[i,], lwd=2, col=tab10.new[i])
  points(to, yo[i,], cex=2, pch=4, col=tab10.new[i])
  abline(h=off[i], lty=3, col=tab10.new[i])
}
legend("topright", legend=c("model", "true", "obs"), 
       lty=c(1, 2, NA), pch=c(NA, NA, 4), ncol=3)
