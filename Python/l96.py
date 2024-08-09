import numpy as np
import rk4

def nl(t, x, F=8):
    return (np.roll(x, -1) - np.roll(x, 2)) * np.roll(x, 1) - x + F

def tl(t, x, dx):
    return (np.roll(dx, -1) - np.roll(dx, 2)) * np.roll(x, 1) + (np.roll(x, -1) - np.roll(x, 2)) * np.roll(dx, 1) - dx

def ad(t, x, dxa):
    return -np.roll(x, -1) * np.roll(dxa, -2) + (np.roll(x, -2) - np.roll(x, 1)) * np.roll(dxa, -1) + np.roll(x, 2) * np.roll(dxa, 1) - dxa

def fom(x, nstep, dt, F):
    t = 0
    for i in range(nstep):
        x = rk4.rk4(nl, t, x, dt, F)
        t += dt
    return x

def tlm(x, dx, nstep, dt):
    t = 0
    for i in range(nstep):
        dx = rk4.tl(nl, tl, t, x, dx, dt)
        t += dt
    return dx

def adm(x, xa, nstep, dt):
    t = nstep
    for i in range(nstep, 0, -1):
        xa = rk4.ad(nl, ad, t, x, xa, dt)
        t -= dt
    return xa

def test(a=1e-5):
    ns = 40
    F = 8
    dt = 0.05
    nstep_spinup = 1000
    seed = 514
    rng = np.random.default_rng(seed)
    
    x0 = fom(rng.normal(size=ns), nstep_spinup, dt, F)
#    x0 = np.ones(ns)
    x = fom(x0, 1, dt, F)
    
    dx0 = a * rng.normal(size=ns)
#    dx0 = a * np.ones(ns)
    x1 = fom(x0 + dx0, 1, dt, F)
    dx = tlm(x0, dx0, 1, dt)
    e = np.sqrt(np.sum((x1 - x - dx)**2))
    print(f"TLM: a={a} l2={e}")
    
    xa = adm(x0, dx, 1, dt)
    tdxdx = np.dot(dx.T, dx)
    dx0xa = np.dot(dx0.T, xa)
    print(f"ADJ: dt^t dx - t(dx0) xa = {tdxdx} - {dx0xa} = {tdxdx - dx0xa}")

# Run the test function
if __name__ == "__main__":
    test()

