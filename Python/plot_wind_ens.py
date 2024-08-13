import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from config_wind import

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
plt.show()

