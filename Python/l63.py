import numpy as np
import ode

def l63(t, w, p, r, b):
  return np.array([
            -p * w[0] + p * w[1],
    (r - w[2]) * w[0] -     w[1],
                     w[0] * w[1] - b * w[2]    
                     ])

def tl(t, wb, wt, p, r, b):
  return np.array([
             -p * wt[0] +     p * wt[1],
    (r - wb[2]) * wt[0] -         wt[1] - wb[0] * wt[2],
          wb[2] * wt[0] + wb[1] * wt[1] -     b * wt[2]
                     ])

def ad(t, wb, dwa, p, r, b):
  return np.array([
    -p * dwa[0] + (r - wb[2]) * dwa[1] + wb[1] * dwa[2],
     p * dwa[0] -               dwa[1] + wb[0] * dwa[2],
                       -wb[0] * dwa[1] -     b * dwa[2]
                     ])

def test(a=1e-5):
    p = 10
    r = 32
    b = 8 / 3
    dt = 0.01
  
    x0 = np.ones(3)
    x = ode.fom(l63, x0, 1, dt, p, r, b)
  
    dx0 = a * x0
    x1 = ode.fom(l63, x0 + dx0, 1, dt, p, r, b)
    dx = ode.tlm(l63, tl, x0, dx0, 1, dt, p, r, b)
    e = np.sqrt(((x1 - x - dx)**2).sum())
    print(f"TLM: a={a} l2={e}")
  
    xa = ode.adm(l63, ad, x0, dx, 1, dt, p, r, b)
    tdxdx = np.dot(dx, dx)
    dx0xa = np.dot(dx0, xa)
    print(f"ADJ: dt^t dx - t(dx0) xa ={tdxdx}-{dx0xa}={tdxdx - dx0xa}")

# Run the test function
if __name__ == "__main__":
    test()
