l96 <- function(t, x, F=8) {
  n <- length(x)
  (x[c(2:n, 1)] - x[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] - x + F
}

l96.tl <- function(t, x, dx){
  n <- length(x)
  (dx[c(2:n, 1)] - dx[c(n-1, n, 1:(n-2))]) * x[c(n, 1:(n-1))] +
    (x[c(2:n, 1)]- x[c(n-1, n, 1:(n-2))]) * dx[c(n, 1:(n-1))] - dx
}

l96.ad <- function(t, xa, dxa){
  n <- length(xa)
  xa <- xa - xa[c(2:n, 1)] * dxa[c(3:n, 1, 2)] + 
    (xa[c(3:n, 1, 2)] - xa[c(n, 1:(n-1))]) * dxa +
    xa[c(n-1, n, 1:(n-2))] * dxa[c(n, 1:(n-1))] - dxa
  dxa = 0.0
  xa
}