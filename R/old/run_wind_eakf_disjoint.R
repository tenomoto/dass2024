set.seed(514)
ne <- 1000
sb <- 2
xb <- c(2, 4)
yo <- 3
so <- 0.3
r <- so^2

calc.speed <- function(u, v) {
  sqrt(u * u + v * v)
}

u <- scale(rnorm(ne, xb[1], sb), scale=FALSE) + xb[1]
v <- scale(rnorm(ne, xb[2], sb), scale=FALSE) + xb[2]
us <- calc.speed(u, v)
xf <- rbind(t(u), t(v))

xa <- xf + eakf.analysis.disjoint(xf, us, yo, r)
xa.mean <- rowMeans(xa)

plot(u, v, xlim = c(-5, 10), ylim = c(-5, 10), asp=1,
     pch=16, cex=0.5, col="gray")
theta <- (0:359) / (2 * pi)
points(xb[1], xb[2], pch=18, cex=1.5, col="black")
points(xa[1,], xa[2,], pch=16, cex=0.5, col="lightblue")
points(xa.mean[1], xa.mean[2], pch=18, cex=1.5, col="blue")

lines((yo + so) * cos(theta), (yo + so) * sin(theta), type="l")
lines((yo - so) * cos(theta), (yo - so) * sin(theta), type="l")
