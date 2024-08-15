source("l63.R")
source("step.R")
source("enkf.R")
source("eakf.R")
source("config_l63.R")

seed <- 514
set.seed(seed)

xt <- matrix(rep(0, ns * nt), nrow=ns)
xt[, 1] <- xt0
yo <- matrix(rep(0, ns * nobs), nrow=ns)
m = 1
for (k in 2:nt) {
  xt[, k] <- step.fom(l63, xt[, k-1], 1, dt, p, r, b)
  if (k %% obs.int == 0) {
    yo[, m] <- rnorm(ns, xt[, k], sqrt(obs.r))
    m <- m + 1
  }
}

xf <- matrix(rep(0, ns * ne), nrow=ns)
for (j in 1: ne) {
  xf[, j] <- rnorm(ns, xb0, sqrt(mstep..q))
}
m = 1
xm <- matrix(rep(0, ns*nt), nrow=ns)
for (k in 1:nt) {
  if (k %% obs.int == 0) {
    dz <- matrix(rep(0, (ns + 1) * ne), ncol=ne)
    for (i in 1:ns) {
      zf <- rbind(xf, xf[i,])
      if (fil == "enkf") {
        yo.p <- rnorm(ns, yo[i, m], sqrt(obs.r[i]))
        dz <- dz + enkf.analysis(zf, yo.p, obs.r[i])
      } else {
        dz <- dz + eakf.analysis(zf, yo[i, m], obs.r[i])
      }
    }
    xf <- xf + dz[1:ns,]
    m <- m + 1
  }
  xm[, k] <- apply(xf, 1, mean)
  if (k < nt){
    for (j in 1:ne) {
      xf[, j] <- step.fom(l63, xf[, j], 1, dt, p, r, b)
    }    
  }
}

tab10.new <- c('#5778a4', '#e49444', '#d1615d', '#85b6b2', '#6a9f58',
               '#e7ca60', '#a87c9f', '#f1a2a9', '#967662', '#b8b0ac')

title <- paste("L63", fil, " ne=", ne)
off <- c(0, 40, 60)

t <- (1:nt) * dt
to <- ((1:nobs) * obs.int) * dt
xt <- xt + off
xm <- xm + off
yo <- yo + off
plot(1, type = "n", xlab="time", ylab="state", main=title,
     xlim=c(min(t), max(t)), ylim=c(-10, 120))
for (i in 1:ns) {
  lines(t, xt[i,], lwd=2, lty=2, col=tab10.new[i])
  lines(t, xm[i,], lwd=2, col=tab10.new[i])
  points(to, yo[i,], cex=2, pch=4, col=tab10.new[i])
  abline(h=off[i], lty=3, col=tab10.new[i])
}
legend("topright", legend=c("mean", "true", "obs"), 
       lty=c(1, 2, NA), pch=c(NA, NA, 4), ncol=ns)
