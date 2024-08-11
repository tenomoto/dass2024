# random
seed <- 514

# l96
ns <- 40
F <- 8
dt <- 0.05

# run
nt.spinup <- 1000
nt <- 1000

# var
b <- 1
nw <- 4

# ens
ne <- 20
#fil <- "eakf"
#c.loc <- 4
#infl <- 1.01
fil <- "enkf"
c.loc <- 1
infl <- 1.01

# obs
r <- 4

# opt
ni <- 100
a <- 0.1
g.break <- 1e-8

