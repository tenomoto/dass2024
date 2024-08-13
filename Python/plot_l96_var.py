import numpy as np
import matplotlib.pyplot as plt
from config_l96 import nt, r, nw, a

plt.rcParams["font.size"] = 18

l2 = np.fromfile("l96_var.dat")

title = f"Var step={a}"
fig, ax = plt.subplots()
ax.plot(range(1, nt+1, nw), l2, color="k", label="analysis")
ax.axhline(y=np.sqrt(r), lw=2, ls='--', color="k", label="observation")
ax.axhline(y=np.mean(l2), lw=2, ls=':', color="k",
           label=f"mean {np.mean(l2):.3f}")
ax.set_title(title)
ax.set_xlabel("time")
ax.set_ylabel("L2")
ax.set_ylim(0, 3)
ax.legend(loc="upper left")
plt.show()


