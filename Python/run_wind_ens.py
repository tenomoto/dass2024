import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import eakf
import enkf

seed = 514
rng = np.random.default_rng(seed)
ne = 1000
sf = 2
xf = np.array([2, 4])
yo = 3
so = 0.3
r = so**2
fil = "enkf"

def calc_speed(u, v):
  return np.sqrt(u * u + v * v)

u = rng.normal(xf[0], sf, ne)
u += -u.mean() + xf[0]
v = rng.normal(xf[1], sf, ne)
v += -v.mean() + xf[1]
us = calc_speed(u, v)
zf = np.vstack([u, v, us])
if fil=="enkf":
  yo_p = rng.normal(yo, so, ne)
  yo_p += yo - yo_p.mean()
  za = zf + enkf.analysis(zf, yo_p, r)
else:
    za = zf + eakf.analysis(zf, yo, r)
za_mean = za.mean(axis=1)

fig, ax = plt.subplots()
ax.scatter(u, v, c="gray", s=8)
ax.scatter(xf[0], xf[1], c="k")
ax.scatter(za[0,], za[1,], c="lightblue", s=8)
ax.scatter(za_mean[0], za_mean[1], c="b")
ax.set_aspect(1)
ax.add_patch(Circle([0, 0], yo - so, 
                    fill=False, ec="k"))
ax.add_patch(Circle([0, 0], yo + so,
                    fill=False, ec="k"))
ax.set_title(f"{fil} ne={ne}")
plt.show()

