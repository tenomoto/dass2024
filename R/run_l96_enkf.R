source("l96.R")
source("enkf.R")
source("eakf.R")
source("dist.R")
source("eval.R")

set.seed(514)

ns <- 40
ne <- 20
F <- 8
nt.spinup <- 1000
nt <- 1000

r <- 4
dt <- 0.05

#fil <- "eakf"
#a <- 4
#infl <- 1.01
fil <- "enkf"
a <- 2
infl <- 1.03

xt <- l96.fom(rnorm(ns), nt.spinup, dt, F)
xf <- matrix(rep(0, ns * ne), nrow=ns)
for (i in 1:ne) {
  xf[,i] <- l96.fom(rnorm(ns), nt.spinup, dt , F)
}

e1 <- rep(0, nt)
e2 <- rep(0, nt)
st <- rep(0, nt)
for (k in 1:nt) {
  yo <- rnorm(ns, xt, sqrt(r))
  dz <- matrix(rep(0, ne * (nrow(xf) + 1)), ncol=ne)
  for (j in 1:length(yo)) {
    zf <- rbind(xf, xf[j,])
    d <- linear.periodic.dist(1:ns, j)
    loc.inf <- c(ga(d, a), infl)
    if (fil == "enkf") {
      yo.p <- rnorm(ne, yo[j], sqrt(r))
      yo.p <- yo.p - mean(yo.p) + yo[j]
      dz <- dz + enkf.analysis(zf, yo.p, r, loc.inf)
    } else {
      dz <- dz + eakf.analysis(zf, yo, r, loc.inf)
    }
  }
  xa <- xf + dz[1:ns,]
  e <- calc.error(xa, xt)
  e1[k] <- e$e1
  e2[k] <- e$e2
  st[k] <- mean(apply(xa, 1, sd))
  if (k < nt) {
    xt <- l96.fom(xt, 1, dt, F)
    for (i in 1:ne) {
      xf[,i] <- l96.fom(xa[,i], 1, dt, F)   
    }
  }
}

title <- paste(fil, "ne=", ne, "infl=", infl, "loc=", a)
plot(1:nt, e1, lwd=2,
     main=title, xlab="step", ylab="L2", type="l", ylim=c(0, 3))
abline(h=1.0, lt=2)
lines(1:nt, st, lwd=2, col="blue")
lines(1:nt, e2, lwd=2, col="red")
legend("topleft", legend=c("mean", "spread", "member"),
       lty=1, col=c("black", "blue", "red"))
