enkf.analysis <- function(zf, yo, r, loc.inf=1.0){
  k <- nrow(zf)
  ne <- ncol(zf)
  zf.mean <- apply(zf, 1, mean)
  zf.anom <- sweep(zf, 1, zf.mean)
  s <- loc.inf * as.vector(zf.anom %*% zf.anom[k,]) / (ne - 1)
  outer(s, (yo - zf[k,])) / (r + s[k])
}