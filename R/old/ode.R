rk4 <- function(f, t, y, h, ...) {
  k1 <- f(t, y, ...)
  k2 <- f(t + 0.5 * h, y + 0.5 * h * k1, ...)
  k3 <- f(t + 0.5 * h, y + 0.5 * h * k2, ...)
  k4 <- f(t + h, y + h * k3, ...)
  y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
}

rk4.tl <- function(f, df, t, y, dy, h, ...) {
  k1 <- f(t, y, ...)
  k2 <- f(t + 0.5 * h, y + 0.5 * h * k1, ...)
  k3 <- f(t + 0.5 * h, y + 0.5 * h * k2, ...)
  dk1 <- df(t, y, dy, ...)
  dk2 <- df(t + 0.5 * h, y + 0.5 * h * k1, dy + 0.5 * h * dk1, ...)
  dk3 <- df(t + 0.5 * h, y + 0.5 * h * k2, dy + 0.5 * h * dk2, ...)
  dk4 <- df(t + h, y + h * k3, dy + h * dk3, ...)
  dy + h * (dk1 + 2 * dk2 + 2 * dk3 + dk4) / 6
}

rk4.ad <- function(f, fa, t, y, ya, dya, h, ...) {
  k1a = h * ya / 6
  k2a = h * ya / 3
  k3a = k2a
  k4a = k1a
  k1 <- f(t, y, ...)
  k2 <- f(t + 0.5 * h, y + 0.5 * h * k1, ...)
  k3 <- f(t + 0.5 * h, y + 0.5 * h * k2, ...)
  dfa = fa(t + h, y + h * k3, k4a, ...)
  ya = ya + dfa
  k3a = k3a + h * dfa
  k4a = 0
  dfa = fa(t + 0.5 * h, y + 0.5 * h * k2, k3a, ...)
  ya = ya + dfa
  k2a = k2a + 0.5 * h * dfa
  k3a = 0
  dfa = fa(t + 0.5 * h, y + 0.5 * h * k1, k2a, ...)
  ya = ya + dfa
  k1a = k1a + 0.5 * h * dfa
  k2a = 0
  ya = ya + fa(t, y, k1a, ...)
  k1a = 0
  ya
}