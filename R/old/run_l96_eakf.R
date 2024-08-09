source("l96.R")
source("ode.R")
source("dist.R")
source("eakf.R")


run_forward <- function(x, nstep, dt, F) {
  t <- 0
  for (i in 1:nstep-1) {
    x <- rk4(l96, t, x, dt, F)
    t <- t + dt
  }
  x
}

calc.error <- function(x, xt) {
  dx = x - xt
  e1 <- sqrt(mean(colMeans(dx)^2))
  e2 <- mean(sqrt(rowMeans(dx^2)))
  r <- e1 / e2
  n <- length(xt)
  r.expected <- sqrt(0.5 * (n + 1) / n)
  r.normalized <- r / r.expected
  e <- list(e1, e2, r, r.expected, r.normalized)
  names(e) <- c("e1", "e2", "r", "r.expected", "r.normalized")
  e
}

set.seed(514)

ns <- 40
ne <- 10
F <- 8
nstep.spinup <- 1000
nstep.da <- 1000
step.da.0 <- 200

r <- 4
dt <- 0.05
c.half <- 0.2
a <- 4
infl <- 1.03

xt <- run_forward(rnorm(ns), nstep.spinup, dt, F)
xf <- matrix(rep(0, ns * ne), nrow=ns)
for (i in 1:ne) {
  xf[,i] <- run_forward(rnorm(ns), nstep.spinup, dt , F)
}

xt1.da <- rep(0, nstep.da)
xa1.da <- matrix(rep(0, ne*nstep.da), ncol=ne)
e1.da <- rep(0, nstep.da)
e2.da <- rep(0, nstep.da)
st.da <- rep(0, nstep.da)
for (k in 1:(step.da.0+nstep.da)) {
  if (k > step.da.0) {
    kk <- k - step.da.0
    xt1.da[kk] <- xt[1]
    xa1.da[kk,] <- xa[1,]
    e <- calc.error(xa, xt)
    e1.da[kk] <- e$e1
    e2.da[kk] <- e$e2
    st.da[kk] <- mean(apply(xa, 1, sd))
  }
  yo <- xt + rnorm(ns, 0, sqrt(r))
  dz <- matrix(rep(0, ne * (nrow(xf) + 1)), ncol=ne)
  for (j in 1:length(yo)) {
    zf <- rbind(xf, xf[j,])
    loc.inf <- c(sqrt(2 * pi) * a * dnorm(linear.periodic.dist(1:ns, j), sd=a), infl^2)
#    loc.inf <- c(ga(linear.periodic.dist(1:ns, j), a), infl^2)
    dz <- dz + eakf.analysis(zf, yo[j], r, loc.inf)
  }
  xa <- xf + dz[1:ns,]
  xt <- run_forward(xt, 1, dt, F)
  for (i in 1:ne) {
    xf[,i] <- run_forward(xa[,i], 1, dt, F)   
  }
}

#matplot(step.da.0:nstep.da, xa1.da, type="l")
#lines(step.da.0:nstep.da, xt1.da, col="blue")
#plot(1:ns, xt, type="l")
#points(1:ns, yo)
#lines(1:ns, rowMeans(xf), col="blue")
#lines(1:ns, rowMeans(xa), col="red")
plot((step.da.0+1):(step.da.0+nstep.da), e1.da, type="l")