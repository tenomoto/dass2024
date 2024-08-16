import numpy as np
import l96
import step
from config_l96 import ns, F, nt_spinup, nt, dt, r, ne

seed = 514
rng = np.random.default_rng(seed)

xt = step.fom(l96.nl, rng.normal(0, 1, ns), nt_spinup, dt, F)
xf = np.zeros(ns)
with open("xf_l96.dat", "wb") as f:
    for j in range(ne):
        xf = step.fom(l96.nl, rng.normal(0, 1, ns), nt_spinup, dt, F)
        xf.tofile(f)

f_true = open("xt_l96.dat", "wb")
f_obs = open("yo_l96.dat", "wb")

for k in range(nt):
    yo = xt + rng.normal(0, np.sqrt(r), ns)
    yo.tofile(f_obs)
    xt.tofile(f_true)
    if k < nt:
        xt = step.fom(l96.nl, xt, 1, dt, F)

