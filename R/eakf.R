eakf.analysis <- function(zf, yo, r, loc.inf=1.0){
  k <- nrow(zf)
  ne <- ncol(zf)
  zf.mean <- apply(zf, 1, mean)
  zf.anom <- sweep(zf, 1, zf.mean)
  s <- loc.inf * as.vector(zf.anom %*% zf.anom[k,]) / (ne - 1)
  alpha <- sqrt(r / (r + s[k]))
  dz.mean <- s / (r + s[k]) * (yo - zf.mean[k])
  dz.mean + (alpha - 1) / s[k] * outer(s, zf.anom[k,])
}