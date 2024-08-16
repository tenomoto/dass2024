import numpy as np
import matplotlib.pyplot as plt
import l63
import step
from config_l63 import xt0, xb0, ns, p, r, b, dt, \
        nt, model_q, obs_int, obs_r, nobs, bgd_b, ni, a 
import sys


seed = 514
rng = np.random.default_rng(seed)

def calc_cost(dx, b, dy, r):
    cost = 0.5 * np.dot(dx, dx) / b
    for j in range(dy.shape[1]):
        cost += 0.5 * np.dot(dy[:, j], dy[:, j] / r)
    return cost

xt = np.zeros([ns, nt])
xt[:, 0] = xt0
yo = np.zeros([ns, nobs])
m = 0
for k in range(1,nt):
    xt[:, k] = step.fom(l63.l63, xt[:, k-1], 1, dt, p, r, b)
    if (k+1) % obs_int == 0:
        yo[:, m] = rng.normal(xt[:, k], np.sqrt(obs_r), ns)
        m += 1

l2 = np.zeros(ni)
cost = np.zeros(ni)
xb = np.zeros([ns, nt])
xb[:, 0] = xb0
for i in range(ni):
    m = 0
    for j in range(1, nt):
        xb[:, j] = step.fom(l63.l63, xb[:, j-1], 1, dt, p, r, b)    
    d = xb[:, (np.arange(nobs) + 1) * obs_int - 1] - yo
    ad = np.zeros(ns)
    m = nobs - 1
    for j in range(nt-1, -1, -1):
        ad = step.adm(l63.l63, l63.ad, xb[:, j], ad, 1, dt, p, r, b)
        if (j+1) % obs_int == 0:
            ad += d[:, m] / obs_r
            m -= 1
    xb[:, 0] = xb[:, 0] - a * ad
    l2[i] = np.sqrt(((xb[:, 0] - xt0)**2).mean())
    cost[i]  = calc_cost(xb[:, 0] - xb0, bgd_b, d, obs_r)
    print(f"i={i} L2={l2[i]:.3e} J={cost[i]:.3e}")

tab = ["tab:blue", "tab:orange", "tab:red"]

title = f"L63 Var a={a:.2e} ni={ni} L2={l2[-1]:.3e}"
off = np.array([0, 40, 60])

t = np.arange(1, nt+1) * dt
to = (np.arange(1,nobs+1) * obs_int) * dt
xt = xt + off[:, None]
xb = xb + off[:, None]
yo = yo + off[:, None]
fig, ax = plt.subplots()
for i in range(ns):
    l_xt = ax.plot(t, xt[i,], lw=2, ls="--", c=tab[i])
    l_xb = ax.plot(t, xb[i,], lw=2, c=tab[i])
    l_yo = ax.scatter(to, yo[i,], s=100, marker="x", c=tab[i])
    ax.axhline(off[i], ls=":", c=tab[i])
ax.set_title(title)
ax.set_xlabel("time")
ax.set_ylabel("state")
ax.set_xlim([t[0], t[-1]])
ax.set_ylim([-20, 150])
plt.show()
