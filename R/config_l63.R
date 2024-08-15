# l63
ns <- 3
x0 <- c(1, 3, 5)
x1 <- c(1.1, -3.3, 5.5)
p <- 10
r <- 32
b <- 8 / 3
dt <- 0.01
nt <- 200
model.q <- x0 * 0.1
obs.int <- 60
obs.r <- x0 * 0.1 
nobs <- nt %/% obs.int

# var
b <- 1
ni <- 100
a <- 5e-4

# ens
ne <- 10
fil <- "eakf"
