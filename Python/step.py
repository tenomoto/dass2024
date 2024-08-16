import rk4


def fom(f, x, nstep, dt, *args):
    t = 0
    for i in range(nstep):
        x = rk4.rk4(f, t, x, dt, *args)
        t += dt
    return x

def tlm(f, df, x, dx, nstep, dt, *args):
    t = 0
    for i in range(nstep):
        dx = rk4.tl(f, df, t, x, dx, dt, *args)
        t += dt
    return dx

def adm(f, fa, x, xa, nstep, dt, *args):
    t = (nstep - 1) * dt
    for i in range(nstep-1, -1, -1):
        xa = rk4.ad(f, fa, t, x, xa, dt, *args)
        t -= dt
    return xa
