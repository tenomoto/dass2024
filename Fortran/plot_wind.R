source("config.R")

title <- paste(fil, "ne=", ne)
fname <- paste0("wind_", fil, ".png")
w <- 1024
h <-  800

con <- file("zf.dat", "rb")
zf <- matrix(readBin(con, "numeric", k * ne), nrow=k)
close(con)
con <- file("za.dat", "rb")
za <- matrix(readBin(con, "numeric", k * ne), nrow=k)
close(con)

png(fname, w, h)
plot(zf[1,], zf[2,], xlim = c(-5, 10), 
     xlab = "u", ylab = "v", main = title,
     ylim = c(-5, 10), asp=1,
     cex.lab=2, cex.axis=2, cex.main=3,
     pch=16, cex=1, col="gray")
theta <- (0:359) / (2 * pi)
points(mean(zf[1,]), mean(zf[2,]), pch=18, cex=3, col="black")
points(za[1,], za[2,], pch=16, cex=1, col="lightblue")
points(mean(za[1,]), mean(za[2,]), pch=18, cex=3, col="blue")

lines((yo + so) * cos(theta), (yo + so) * sin(theta), type="l")
lines((yo - so) * cos(theta), (yo - so) * sin(theta), type="l")
dev.off()
