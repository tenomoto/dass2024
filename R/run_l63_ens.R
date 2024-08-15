source("l63.R")
source("ode.R")
source("enkf.R")
source("eakf.R")
source("config_l63.R")

seed <- 514
set.seed(seed)

xt <- matrix(rep(0, ns * nt), nrow=ns)
xt[, 1] <- x0
yo <- matrix(rep(0, ns * nobs), nrow=ns)
m = 1
for (k in 2:nt) {
  xt[, k] <- ode.fom(l63, xt[, k-1], 1, dt, p, r, b)
  if (k %% obs.int == 0) {
    yo[, m] <- rnorm(ns, xt[, k], sqrt(obs.r))
    m <- m + 1
  }
}

xf <- matrix(rep(0, ns * ne), nrow=ns)
for (j in 1: ne) {
  xf[, j] <- rnorm(ns, x1, sqrt(model.q))
}
m = 1
e1 <- rep(0, nobs)
e2 <- rep(0, nobs)
st <- rep(0, nobs)
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
    e <- calc.error(xf, xt[, k])
    e1[m] <- e$e1
    e2[m] <- e$e2
    st[m] <- mean(apply(xf, 1, sd))
    m <- m + 1
  }
  xm[, k] <- apply(xf, 1, mean)
  if (k < nt){
    for (j in 1:ne) {
      xf[, j] <- ode.fom(l63, xf[, j], 1, dt, p, r, b)
    }    
  }
}

title <- paste("L63", fil)
off <- c(0, 40, 60)
t <- (1:nt) * dt
to <- ((1:nobs) * obs.int) * dt
xt <- xt + off
xm <- xm + off
yo <- yo + off
plot(1, type = "n", xlab="time", ylab="x,y,z", main=title,
  xlim=c(min(t), max(t)), ylim=c(min(xm), max(xm)+20))
for (i in 1:ns) {
  lines(t, xt[i,], lwd=2, col="gray")
  lines(t, xm[i,], lwd=2, col="red")
  points(to, yo[i,], cex=2, pch=4, col="blue")
}
legend("topright", legend=c("mean", "true", "obs"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 4),
       col=c("red", "gray", "blue"), ncol=ns)
