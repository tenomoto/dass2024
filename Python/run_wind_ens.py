import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import eakf

seed = 514
rng = np.random.default_rng(seed)
ne = 1000
sf = 2
xf = np.array([2, 4])
yo = 3
so = 0.3
r = so**2

def calc_speed(u, v):
  return np.sqrt(u * u + v * v)

u = rng.normal(xf[0], sf, ne)
u += -u.mean() + xf[0]
v = rng.normal(xf[1], sf, ne)
v += -v.mean() + xf[1]
us = calc_speed(u, v)
zf = np.vstack([u, v, us])
za = zf + eakf.analysis(zf, yo, r)
za_mean = za.mean(axis=1)

