source("l63.R")
source("step.R")

ns <- 3
x0 <- c(1, 3, 5)
x0.p <- c(1.1, 3.3, 5.5)
p <- 10
r <- 48
b <- 8 / 3
dt <- 0.01
nt <- 500

x.rk4 <- matrix(rep(0, ns * nt), nrow=ns)
x.rk4[, 1] <- x0
x.rk4p <- matrix(rep(0, ns * nt), nrow=ns)
x.rk4p[, 1] <- x0.p
x.eul <- matrix(rep(0, ns * nt), nrow=ns)
x.eul[, 1] <- x0

for (j in 2:nt) {
  t <- (j - 1) * dt
  x.rk4[, j] <- step.fom(l63, x.rk4[, j-1], 1, dt, p, r, b)
  x.rk4p[, j] <- step.fom(l63, x.rk4p[, j-1], 1, dt, p, r, b) 
  x.eul[, j] <- x.eul[, j-1] + l63(t, x.eul[, j-1], p, r, b) * dt
}

tab10.new <- c('#5778a4', '#e49444', '#d1615d', '#85b6b2', '#6a9f58',
               '#e7ca60', '#a87c9f', '#f1a2a9', '#967662', '#b8b0ac')

title <- paste("L63 free run")
off <- c(0, 40, 60)
t <- (1:nt) * dt
x.rk4 <- x.rk4 + off
x.rk4p <- x.rk4p + off
x.eul <- x.eul + off
plot(1, type = "n", xlab="time", ylab="state", main=title,
     xlim=c(min(t), max(t)), ylim=c(-20, 150))
for (i in 1:ns) {
  lines(t, x.rk4[i,], lwd=2, lty=1, col=tab10.new[i])
  lines(t, x.rk4p[i,], lwd=2, lty=2, col=tab10.new[i])
  lines(t, x.eul[i,], lwd=2, lty=3, col=tab10.new[i])
  abline(h=off[i], lty=3, col=tab10.new[i])
}
legend("topright", legend=c("RK4", "RK4p", "Euler"), lty=c(1, 2, 3), ncol=3)
