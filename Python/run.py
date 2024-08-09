import numpy as np
import l96

seed = 514
rng = np.random.default_rng(seed)

ns = 40
F = 8

nw = 4
nt_spinup = 1000
nt = 1000
nc = nt // nw
ni = 100

r = 4
b = 1
dt = 0.05
a = 0.1
g_break = 1e-8

def calc_cost(dx, b, dy, r):
    cost = 0.5 * np.dot(dx, dx) / b
    for i in range(dy.shape[1]):
        cost += 0.5 * np.dot(dy[:, j], dy[:, j]) / 4
    return cost

xt = l96.fom(rng.normal(size=ns), nt_spinup, dt, F)
x0 = l96.fom(rng.normal(size=ns), nt_spinup, dt, F)

yo = np.zeros([ns, nw])
xb = np.zeros([ns, nw])
d = np.zeros([ns, nw])
e = np.zeros(nc)

print(xt)
print(x0)
print(np.sqrt(np.mean((xt-x0)**2)))

for k in range(nc):
    print(f"cycle={k+1}")
    for j in range(nw):
        yo[:, j] = rng.normal(xt, np.sqrt(r), ns)
        xt = l96.fom(xt, 1, dt, F)
    xb[:, 0] = x0[:]
    for in in range(ni):
        for j in range(1, nw):
            xb[:, j] = l96.adm(xb[:, j], 1, dt, F)
        d = xb - yo
        ad = np.zeros(ns)
        for j in range(nw, -1, -1):
            ad = l96.adm(xb[:, j], ad, 1, dt) + d[:, j] / r
        xb[:, 0] -= a * ad
        cost = calc_cost(xb[:, 0] - x0, b, d, r)
        if i==1: cost.old = cost
        e[k] = np.sqrt(np.mean(xb[:, 0] - xt)**2

