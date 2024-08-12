source("enkf.R")
source("eakf.R")

set.seed(514)
ne <- 1000
sf <- 2
xf <- c(2, 4)
yo <- 3
so <- 0.3
r <- so^2
fil <- "eakf"

calc.speed <- function(u, v) {
  sqrt(u * u + v * v)
}

u <- scale(rnorm(ne, xf[1], sf), scale=FALSE) + xf[1]
v <- scale(rnorm(ne, xf[2], sf), scale=FALSE) + xf[2]
us <- calc.speed(u, v)
zf <- rbind(t(u), t(v), t(us))
if (fil=="enkf") {
  yo.p <- rnorm(ne, yo, so)
  yo.p <- yo.p - mean(yo.p) + yo
  za <- zf + enkf.analysis(zf, yo.p, r)
} else {
  za <- zf + eakf.analysis(zf, yo, r)
}
za.mean <- apply(za, 1, mean)

plot(u, v, xlim = c(-5, 10), ylim = c(-5, 10), asp=1,
     pch=16, cex=0.5, col="gray")
theta <- (0:359) / (2 * pi)
points(xf[1], xf[2], pch=18, cex=1.5, col="black")
points(za[1,], za[2,], pch=16, cex=0.5, col="lightblue")
points(za.mean[1], za.mean[2], pch=18, cex=1.5, col="blue")

lines((yo + so) * cos(theta), (yo + so) * sin(theta), type="l")
lines((yo - so) * cos(theta), (yo - so) * sin(theta), type="l")
