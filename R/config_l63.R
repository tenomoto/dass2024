# l63
ns <- 3
xt0 <- c(1, 3, 5)
#xb0 <- c(1.1, 3.3, 5.5)
xb0 <- c(1, -1, 1)
p <- 10
r <- 48
b <- 8 / 3
dt <- 0.01
nt <- 200
model.q <- 0.1
obs.int <- 60
obs.r <- c(0.1, 0.3, 0.5)^2
nobs <- nt %/% obs.int

# var
bgd.b <- 1.0
ni <- 10
a <- 5e-4

# ens
ne <- 10
fil <- "eakf"
