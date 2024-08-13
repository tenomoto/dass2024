import numpy as np
import l96
import enkf
import eakf
import dist
from eval import calc_error
from config_l96 import ns, F, nt_spinup, nt, dt, r, \
        ne, fil, c_loc, infl, xt_fname, yo_fname, xf_fname

seed = 516
rng = np.random.default_rng(seed)

f_xt = open(xt_fname, "rb")
f_yo = open(yo_fname, "rb")
xf = np.zeros([ns, ne])
with open("xf_l96.dat", "rb") as f:
    for j in range(ne):
        xf[:, j] = np.fromfile(f, count=ns)

e1 = np.zeros(nt)
e2 = np.zeros(nt)
st = np.zeros(nt)
for k in range(nt):
    xt = np.fromfile(f_xt, count=ns)
    yo = np.fromfile(f_yo, count=ns)
    dz = np.zeros([ns+1, ne])
    for i in range(ns):
        zf = np.vstack([xf, xf[i,]])
        d = dist.linear_periodic(np.arange(ns), i)
        loc_inf = np.append(dist.ga(d, c_loc), infl)
        if fil == "enkf":
            yo_p = rng.normal(yo[i], np.sqrt(r), ne)
            yo_p += yo[i] - yo_p.mean()
            dz = dz + enkf.analysis(zf, yo_p, r, loc_inf)
        else:
            dz = dz + eakf.analysis(zf, yo[i], r, loc_inf)
    xa = xf + dz[0:-1, :]
    e = calc_error(xa, xt)
    e1[k] = e[0]
    e2[k] = e[1]
    st[k] = xa.std(axis=1).mean()
    if k < nt - 1:
        for i in range(ne):
            xf[:, i] = l96.fom(xa[:, i], 1, dt, F)
f_xt.close()
f_yo.close()

with open(f"l96_{fil}.dat", "wb") as f_out:
    e1.tofile(f_out)
    e2.tofile(f_out)
    st.tofile(f_out)
