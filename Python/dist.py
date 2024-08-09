import numpy as np


def lorenc(i, j, n=40) :
    return np.abs(n/np.pi * np.sin(np.pi * (i - j) / n))

def linear_periodic(i, j, n=40):
    return np.minimum(np.abs(i - j), n - np.abs(i - j))

def gc_5th(z):
    r = np.abs(z)
    return np.where(r <= 1,
        -0.25 * r**5 + 0.5 * r**4 + 5/8 * r**3 - 5/3 * r**2 + 1,
        np.where (r <= 2,
            1/12 * r**5 - 0.5 * r**4 + 5/8 * r**3 + 5/3 * r**2 - 5 * r + 4 -2/3/r,
            0))

def ga(r, a):
    return np.exp(-0.5 * (r / a)**2)
