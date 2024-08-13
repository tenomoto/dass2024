import numpy as np


def calc_error(x, xt):
    dx = x - xt[:, None]
    e1 = np.sqrt((dx.mean(axis=1)**2).mean())
    e2 = (np.sqrt((dx**2).mean(axis=0))).mean()
    r = e1 / e2
    n = xt.size
    r_expected = np.sqrt(0.5 * (n + 1) / n)
    r_normalized = r / r_expected
    return np.array([e1, e2, r, r_expected, r_normalized])
