import numpy as np
import l96
from config_l96 import ns, F, dt, r, nt, \
        b, nw, ni, a, g_break, xt_fname, yo_fname, xf_fname

nc = nt // nw
def calc_cost(dx, b, dy, r):
    cost = 0.5 * np.dot(dx, dx) / b
    for j in range(dy.shape[1]):
        cost += 0.5 * np.dot(dy[:, j], dy[:, j]) / r
    return cost

f_xt = open(xt_fname, "rb")
f_yo = open(yo_fname, "rb")
x0 = np.fromfile(xf_fname, count=ns)

yo = np.zeros((ns, nw))
xb = np.zeros((ns, nw))
d = np.zeros((ns, nw))
l2 = np.zeros(nc)
ad = np.zeros(ns)

for k in range(nc):
    print(f"cycle={k+1}")
    xt0 = np.fromfile(f_xt, count=ns)
    for j in range(1, nw):
        xt = np.fromfile(f_xt, count=ns)
    for j in range(0, nw):
        yo[:, j] = np.fromfile(f_yo, count=ns)
    xb[:, 0] = x0
    for i in range(ni):
        for j in range(1, nw):
            xb[:, j] = l96.fom(xb[:, j-1], 1, dt, F)
        d = xb - yo
        ad = 0
        for j in range(nw-1, -1, -1):
            ad = l96.adm(xb[:, j], ad, 1, dt, F) + d[:, j] / r
        xb[:, 0] -= a * ad
        cost = calc_cost(xb[:, 0] - x0, b, d, r)
        if i == 0:
            cost_old = cost
        l2[k] = np.sqrt(np.mean((xb[:, 0] - xt0)**2))
        gnorm = np.dot(ad, ad)
        print(f"{i+1} J={cost:.2e} Jold={cost_old:.2e} g={gnorm:.2e} l2={l2[k]:.2e}")
        if gnorm < g_break:
            break
        if cost > cost_old:
            xb[:, 0] += a * ad
            l2[k] = np.sqrt(np.mean((xb[:, 0] - xt0)**2))
            break
        cost_old = cost
    x0 = l96.fom(xb[:, 0], nw + 1, dt, F)
f_xt.close()
f_yo.close()

l2.tofile("l96_var.dat")
