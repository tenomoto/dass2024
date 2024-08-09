calc.error <- function(x, xt) {
  dx = sweep(x, 1, xt)
  e1 <- sqrt(mean(apply(dx, 1, mean)^2))
  e2 <- mean(sqrt(apply(dx^2, 2, mean)))
  r <- e1 / e2
  n <- length(xt)
  r.expected <- sqrt(0.5 * (n + 1) / n)
  r.normalized <- r / r.expected
  e <- list(e1, e2, r, r.expected, r.normalized)
  names(e) <- c("e1", "e2", "r", "r.expected", "r.normalized")
  e
}