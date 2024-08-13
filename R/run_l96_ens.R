source("l96.R")
source("enkf.R")
source("eakf.R")
source("dist.R")
source("eval.R")
source("config_l96.R")

seed <- 516
set.seed(seed)

con.true <- file("xt_l96.dat", "rb")
con.obs <- file("yo_l96.dat", "rb")
con.ens <- file("xf_l96.dat", "rb")
xf <- matrix(readBin(con.ens, "numeric", ns * ne), nrow=ns)
close(con.ens)

e1 <- rep(0, nt)
e2 <- rep(0, nt)
st <- rep(0, nt)
for (k in 1:nt) {
  xt <- readBin(con.true, "numeric", ns)
  yo <- readBin(con.obs, "numeric", ns)
  dz <- matrix(rep(0, (ns + 1) * ne), ncol=ne)
  for (i in 1:ns) {
    zf <- rbind(xf, xf[i,])
    d <- linear.periodic.dist(1:ns, i)
    loc.inf <- c(ga(d, c.loc), infl)
    if (fil == "enkf") {
      yo.p <- rnorm(ne, yo[i], sqrt(r))
      yo.p <- yo.p - mean(yo.p) + yo[i]
      dz <- dz + enkf.analysis(zf, yo.p, r, loc.inf)
    } else {
      dz <- dz + eakf.analysis(zf, yo[i], r, loc.inf)
    }
  }
  xa <- xf + dz[1:ns,]
  e <- calc.error(xa, xt)
  e1[k] <- e$e1
  e2[k] <- e$e2
  st[k] <- mean(apply(xa, 1, sd))
  if (k < nt) {
    for (i in 1:ne) {
      xf[,i] <- l96.fom(xa[,i], 1, dt, F)   
    }
  }
}
close(con.true)
close(con.obs)

fname <- paste0("l96_", fil, ".dat")
con.out <- file(fname, "wb")
writeBin(e1, con.out)
writeBin(e2, con.out)
writeBin(st, con.out)
close(con.out)
