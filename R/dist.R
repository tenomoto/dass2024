lorenc.dist <- function(i, j, n=40) {
  abs(n/pi * sin(pi * (i - j) / n))
}

linear.periodic.dist <- function(i, j, n=40) {
  pmin(abs(i - j), n - abs(i - j))
}

gc.5th <- function(z) {
  r <- abs(z)
  ifelse (r <= 1,
          -0.25 * r^5 + 0.5 * r^4 + 5/8 * r^3 - 5/3 * r^2 + 1,
          ifelse (r <= 2,
                  1/12 * r^5 - 0.5 * r^4 + 5/8 * r^3 + 5/3 * r^2 - 5 * r + 4 -2/3/r,
                  0)
  )
}

ga <- function(r, a) {
  sqrt(2 * pi) * a * dnorm(r, sd=a)
}