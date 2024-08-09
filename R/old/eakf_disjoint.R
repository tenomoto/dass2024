eakf.analysis.disjoint <- function(xf, yf, yo, r){
  k <- nrow(xf)
  ne <- ncol(xf)
  xf.anom <- t(scale(t(xf), scale=FALSE))
  xf.mean <- attr(xf.anom, "scaled:center")
  yf.mean <- mean(yf)
  
  yf.anom <- as.vector(yf - yf.mean)
  sxy <- as.vector(xf.anom %*% yf.anom) / (ne - 1)
  syy <- as.numeric(t(yf.anom) %*% yf.anom) / (ne - 1)
  alpha <- sqrt(r / (r + syy))
  dz.mean <- sxy / (r + syy) * (yo - yf.mean)
  dz <-  dz.mean + (alpha - 1) / syy * outer(sxy, yf.anom)
}