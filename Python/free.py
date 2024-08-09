import numpy as np
import matplotlib.pyplot as plt
import l96
import rk4

nj = 40
nstep = 1001
x_hist = np.zeros((nj, nstep))
x = np.random.normal(size=nj)
F = 3.85
dt = 0.05
t = 0.0

for i in range(nstep - 1):
    x = rk4.rk4(l96.nl, t, x, dt, F)
    x_hist[:, i] = x
    t += dt

t = np.linspace(0, nstep * dt, nstep)

fig, ax = plt.subplots()
cnt = ax.contourf(range(1, nj + 1), t, x_hist.T,
                  levels=11, cmap="YlOrRd")
ax.invert_yaxis()
ax.set_xlabel('j')
ax.set_ylabel('time')
fig.colorbar(cnt)
plt.show()
