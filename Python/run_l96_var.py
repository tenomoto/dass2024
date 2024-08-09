import numpy as np
import matplotlib.pyplot as plt
import l96
import sys

def calc_cost(dx, b, dy, r):
    cost = 0.5 * np.dot(dx, dx) / b
    for j in range(dy.shape[1]):
        cost += 0.5 * np.dot(dy[:, j], dy[:, j]) / r
    return cost

seed = 514
rng = np.random.default_rng(seed)

ns = 40
F = 8

nw = 4
nt_spinup = 1000
nt = 1000
nc = nt // nw
ni = 1000

r = 4
b = 1
dt = 0.05
a = 0.1
g_break = 1e-8

# Initial conditions
xt = l96.fom(rng.normal(size=ns), nt_spinup, dt, F)
x0 = l96.fom(rng.normal(size=ns), nt_spinup, dt, F)

yo = np.zeros((ns, nw))
xb = np.zeros((ns, nw))
d = np.zeros((ns, nw))
e = np.zeros(nc)

for k in range(nc):
    print(f"cycle={k+1}")
    xt0 = xt.copy()
    for j in range(nw):
        yo[:, j] = rng.normal(loc=xt, scale=np.sqrt(r), size=ns)
        xt = l96.fom(xt, 1, dt, F)
    
    xb[:, 0] = x0
    for i in range(ni):
        for j in range(1, nw):
            xb[:, j] = l96.fom(xb[:, j-1], 1, dt, F)
        
        d = xb - yo
        ad = np.zeros(ns)
        for j in range(nw-1, -1, -1):
            ad = l96.adm(xb[:, j], ad, 1, dt) + d[:, j] / r
        
        xb[:, 0] -= a * ad
        cost = calc_cost(xb[:, 0] - x0, b, d, r)
        if i == 0:
            cost_old = cost
        
        e[k] = np.sqrt(np.mean((xb[:, 0] - xt0)**2))
        gnorm = np.dot(ad, ad)
        print(f"{i+1} J={cost:.2e} Jold={cost_old:.2e} g={gnorm:.2e} l2={e[k]:.2e}")
        
        if gnorm < g_break:
            break
        if cost > cost_old:
            xb[:, 0] += a * ad
            break
        cost_old = cost
    
    x0 = l96.fom(xb[:, 0], nw + 1, dt, F)
#fig, ax = plt.subplots()
#ax.plot(xt)
#ax.plot(x0)
#ax.plot(xb[:,1])
#ax.plot(xb[:,1] - x0)
#plt.show()


title = f"4DVar step={a}"
fig, ax = plt.subplots()
ax.plot(range(1, nt+1, nw), e, color="k", label="analysis")
ax.axhline(y=1, linestyle='--', color="k", label="observation")
ax.axhline(y=np.mean(e), linestyle=':', color="k", label=f"mean {np.mean(e):.3f}")
ax.set_title(title)
ax.set_xlabel("time")
ax.set_ylabel("error")
ax.set_ylim(0, 3)
ax.legend(loc="upper left")
plt.show()

