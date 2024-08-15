source("l96.R")
source("step.R")
source("config_l96.R")

seed <- 514
set.seed(seed)

xt <- step.fom(l96, rnorm(ns), nt.spinup, dt, F)
con.ens <- file(xf_fname, "wb")
for (j in 1:ne) {
  xf <- step.fom(l96, rnorm(ns), nt.spinup, dt , F)
  writeBin(xf, con.ens)
}
close(con.ens)

con.true <- file(xt_fname, "wb")
con.obs <- file(yo_fname, "wb")
for (k in 1:nt) {
  yo <- rnorm(ns, xt, sqrt(r))
  writeBin(yo, con.obs)
  writeBin(xt, con.true)
  if (k < nt) {
    xt <- step.fom(l96, xt, 1, dt, F)
  }
}
close(con.true)
close(con.obs)
