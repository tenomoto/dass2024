source("ode.R")

l63 <- function(t, w, p, r, b){
  c(
            -p * w[1] + p * w[2],
    (r - w[3]) * w[1] -     w[2],
                     w[1] * w[2] - b * w[3]    
  )
}

l63.tl <- function(t, wb, wt, p, r, b){
  c(
             -p * wt[1] +     p * wt[2],
    (r - wb[3]) * wt[1] -         wt[2] - wb[1] * wt[3],
          wb[2] * wt[1] + wb[1] * wt[2] -     b * wt[3]
  )
}

l63.ad <- function(t, wb, dwa, p, r, b){
  c(
    -p * dwa[1] + (r - wb[3]) * dwa[2] + wb[2] * dwa[3],
     p * dwa[1] -               dwa[2] + wb[1] * dwa[3],
                       -wb[1] * dwa[2] -     b * dwa[3]
  )
}

l63.test <- function(a=1e-5) {
  p <- 10
  r <- 32
  b <- 8 / 3
  dt <- 0.01
  
  x0 <- c(1, 1, 1)
  x <- ode.fom(l63, x0, 1, dt, p, r, b)
  
  dx0 <- a * x0
  x1 <- ode.fom(l63, x0 + dx0, 1, dt, p, r, b)
  dx <- ode.tlm(l63, l63.tl, x0, dx0, 1, dt, p, r, b)
  e <- sqrt(sum((x1 - x - dx)^2))
  cat("TLM:", "a=", a, "l2=", e, "\n")
  
  xa <- ode.adm(l63, l63.ad, x0, dx, 1, dt, p, r, b)
  tdxdx <- t(dx) %*% dx
  dx0xa <- t(dx0) %*% xa
  cat("ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa, "\n")
}
