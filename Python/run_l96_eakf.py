import numpy as np
import matplotlib.pyplot as plt
import l96
import rk4
import eakf
import dist


def calc_error(x, xt):
    dx = x - xt[:, None]
    e1 = np.sqrt((dx.mean(axis=0)**2).mean())
    e2 = np.sqrt((dx**2).mean(axis=1)).mean()
    r = e1 / e2
    n = xt.size
    r_expected = np.sqrt(0.5 * (n + 1) / n)
    r_normalized = r / r_expected
    return np.array([e1, e2, r, r_expected, r_normalized])

seed = 514
rng = np.random.default_rng(seed)

ns = 40
ne = 10
F = 8
nstep_spinup = 1000
nstep_data = 1000
step_da_0 = 200

r = 4
dt = 0.05
c_half = 0.3
infl = 1.03

xt = l96.fom(rng.normal(0, 1, ns), nstep_spinup, dt, F)
xf = np.zeros([ns, ne])
for i in range(ne):
    xf[:,i] = l96.fom(rng.normal(0, 1, ns), nstep_spinup, dt, F)

xt1_da = np.zeros(nstep_data)
xa1_da = np.zeros([nstep_data, ne])
e1_da = np.zeros(nstep_data)
e2_da = np.zeros(nstep_data)
st_da = np.zeros(nstep_data)
for k in range(step_da_0 + nstep_data):
    if k > step_da_0:
        kk = k - step_da_0
        xt1_da[kk] = xt[1]
        xa1_da[kk, :] = xa[1, :]
        e = calc_error(xa, xt)
        e1_da[kk] = e[0]
        e2_da[kk] = e[1]
        st_da[kk] = xa.std(axis=1).mean()
    yo = xt + rng.normal(0, np.sqrt(r), ns)
    dz = np.zeros([ns+1, ne])
    for j in range(ne):
        zf = np.vstack([xf, xf[j,]])
        loc_inf = np.append(
##            dist.gc_5th(dist.linear_periodic(np.arange(ns), j) /(c_half*ns/2)),
            dist.ga(dist.linear_periodic(np.arange(ns), j), 5),
            infl**2)
        dz = dz + eakf.analysis(zf, yo[j], r, loc_inf)
    xa = xf + dz[0:-1]
    xt = l96.fom(xt, 1, dt, F)
    for i in range(ne):
        xf[:, i] = l96.fom(xa[:, i], 1, dt, F)   

fig, ax = plt.subplots()
ax.plot(np.arange(step_da_0+1, step_da_0+nstep_data+1), e1_da)
plt.show()
