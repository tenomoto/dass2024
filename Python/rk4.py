def rk4(f, t, y, h, *args):
    k1 = f(t, y, *args)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1, *args)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2, *args)
    k4 = f(t + h, y + h * k3, *args)
    return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

def tl(f, df, t, y, dy, h, *args):
    k1 = f(t, y, *args)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1, *args)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2, *args)
    dk1 = df(t, y, dy, *args)
    dk2 = df(t + 0.5 * h, y + 0.5 * h * k1, dy + 0.5 * h * dk1, *args)
    dk3 = df(t + 0.5 * h, y + 0.5 * h * k2, dy + 0.5 * h * dk2, *args)
    dk4 = df(t + h, y + h * k3, dy + h * dk3, *args)
    return dy + h * (dk1 + 2 * dk2 + 2 * dk3 + dk4) / 6

def ad(f, fa, t, y, ya, h, *args):
    k1 = f(t, y, *args)
    y2 = y + 0.5 * h * k1
    k2 = f(t + 0.5 * h, y2, *args)
    y3 = y + 0.5 * h * k2
    k3 = f(t + 0.5 * h, y3, *args)
    y4 = y + h * k3

    k1a = h * ya / 6
    k2a = h * ya / 3
    k3a = k2a
    k4a = k1a
    
    dfa = fa(t + h, y4, k4a, *args)
    ya = ya + dfa
    k3a = k3a + h * dfa
    k4a = 0
    dfa = fa(t + 0.5 * h, y3, k3a, *args)
    ya = ya + dfa
    k2a = k2a + 0.5 * h * dfa
    k3a = 0
    dfa = fa(t + 0.5 * h, y2, k2a, *args)
    ya = ya + dfa
    k1a = k1a + 0.5 * h * dfa
    k2a = 0
    ya = ya + fa(t, y, k1a, *args)
    k1a = 0
    return ya

