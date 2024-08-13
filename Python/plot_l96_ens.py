import numpy as np
import matplotlib.pyplot as plt
from config_l96 import nt, r, ne, fil, c_loc, infl


with open(f"l96_{fil}.dat", "rb") as f:
    e1 = np.fromfile(f, count=nt)
    e2 = np.fromfile(f, count=nt)
    st = np.fromfile(f, count=nt)

plt.rcParams["font.size"] = 18

title = f"{fil} ne={ne} infl={infl} loc={c_loc}"
fig, ax = plt.subplots()
ax.plot(np.arange(nt) + 1, e1, lw=2, c="black")
ax.plot(np.arange(nt) + 1, e2, lw=2, c="red")
ax.plot(np.arange(nt) + 1, st, lw=2, c="blue")
ax.set_title(title)
ax.set_xlabel("time")
ax.set_ylabel("L2")
ax.set_ylim([0, 3])
ax.axhline(np.sqrt(r), lw=2, ls="--", c="black")
plt.show()
