source("l96.R")
source("config_l96.R")

seed <- 514
set.seed(seed)

xt <- l96.fom(rnorm(ns), nt.spinup, dt, F)
con.ens <- file("xf_l96.dat", "wb")
for (i in 1:ne) {
  xf <- l96.fom(rnorm(ns), nt.spinup, dt , F)
  writeBin(xf, con.ens)
}
close(con.ens)

con.true <- file("xt_l96.dat", "wb")
con.obs <- file("yo_l96.dat", "wb")
for (k in 1:nt) {
  yo <- rnorm(ns, xt, sqrt(r))
  writeBin(yo, con.obs)
  writeBin(xt, con.true)
  if (k < nt) {
    xt <- l96.fom(xt, 1, dt, F)
  }
}
close(con.true)
close(con.obs)
