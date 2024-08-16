import numpy as np
import step

def nl(t, x, F=8):
    return (np.roll(x, -1) - np.roll(x, 2)) * np.roll(x, 1) - x + F

def tl(t, x, dx, *args):
    return (np.roll(dx, -1) - np.roll(dx, 2)) * np.roll(x, 1) + (np.roll(x, -1) - np.roll(x, 2)) * np.roll(dx, 1) - dx

def ad(t, x, dxa, *args):
    return -np.roll(x, -1) * np.roll(dxa, -2) + (np.roll(x, -2) - np.roll(x, 1)) * np.roll(dxa, -1) + np.roll(x, 2) * np.roll(dxa, 1) - dxa

def test(a=1e-5):
    ns = 40
    F = 8
    dt = 0.05
    nstep_spinup = 1000
    seed = 514
    rng = np.random.default_rng(seed)
    
    x0 = step.fom(nl, rng.normal(size=ns), nstep_spinup, dt, F)
#    x0 = np.ones(ns)
    x = step.fom(nl, x0, 1, dt, F)
    
    dx0 = a * rng.normal(size=ns)
#    dx0 = a * np.ones(ns)
    x1 = step.fom(nl, x0 + dx0, 1, dt, F)
    dx = step.tlm(nl, tl, x0, dx0, 1, dt, F)
    e = np.sqrt(((x1 - x - dx)**2).mean())
    print(f"TLM: a={a} l2={e}")
    
    xa = step.adm(nl, ad, x0, dx, 1, dt, F)
    tdxdx = np.dot(dx.T, dx)
    dx0xa = np.dot(dx0.T, xa)
    print(f"ADJ: dt^t dx - t(dx0) xa = {tdxdx} - {dx0xa} = {tdxdx - dx0xa}")

# Run the test function
if __name__ == "__main__":
    test()

