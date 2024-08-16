import numpy as np
import matplotlib.pyplot as plt
from l63 import l63
import step
import enkf
import eakf
from config_l63 import xt0, xb0, ns, p, r, b, dt, \
        nt, model_q, obs_int, obs_r, nobs, ne, fil


seed = 514
rng = np.random.default_rng(seed)

xt = np.zeros([ns, nt])
xt[:, 0] = xt0
yo = np.zeros([ns, nobs])
m = 0
for k in range(1, nt):
    xt[:, k] = step.fom(l63, xt[:, k-1], 1, dt, p, r, b)
    if (k+1) % obs_int == 0:
        yo[:, m] = rng.normal(xt[:, k], np.sqrt(obs_r), ns)
        m += 1

xf = np.zeros([ns, ne])
for i in range(ns):
    xf[i, :] = rng.normal(xb0[i], np.sqrt(model_q), ne)
m = 0
xm = np.zeros([ns, nt])
for k in range(nt):
    if (k + 1) % obs_int == 0:
        dz = np.zeros([ns+1, ne])
        for i in range(ns):
            zf = np.vstack([xf, xf[i,:]])
            if fil == "enkf":
                yo_p = rng.normal(yo[i, m], np.sqrt(obs_r[i]), nw)
                dz += enkf.analysis(zf, yo_p, obs_r[i])
            else:
                dz += eakf.analysis(zf, yo[i, m], obs_r[i])
        xf += dz[0:ns,:]
        m += 1
    xm[:, k] = xf.mean(axis=1)
    if k < nt:
        for j in range(ne):
            xf[:, j] = step.fom(l63, xf[:, j], 1, dt, p, r, b)

tab = ["tab:blue", "tab:orange", "tab:red"]

title = f"L63 {fil} ne={ne}"
off = np.array([0, 40, 60])

t = np.arange(1, nt+1) * dt
to = (np.arange(1,nobs+1) * obs_int) * dt
xt = xt + off[:, None]
xm = xm + off[:, None]
yo = yo + off[:, None]
fig, ax = plt.subplots()
for i in range(ns):
    l_xt = ax.plot(t, xt[i,], lw=2, ls="--", c=tab[i])
    l_xm = ax.plot(t, xm[i,], lw=2, c=tab[i])
    l_yo = ax.scatter(to, yo[i,], s=100, marker="x", c=tab[i])
    ax.axhline(off[i], ls=":", c=tab[i])
ax.set_title(title)
ax.set_xlabel("time")
ax.set_ylabel("state")
ax.set_xlim([t[0], t[-1]])
ax.set_ylim([-20, 150])
plt.show()
