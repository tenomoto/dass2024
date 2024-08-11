source("config_l96.R")

fname <- paste0("l96_", fil, ".dat")
con <- file(fname, "rb")
e1 <- readBin(con, "numeric", nt)
e2 <- readBin(con, "numeric", nt)
st <- readBin(con, "numeric", nt)
close(con)

title <- paste(fil, "ne=", ne, "infl=", infl, "loc=", c.loc)
plot(1:nt, e1, lwd=2,
     main=title, xlab="step", ylab="L2", type="l", ylim=c(0, 3))
abline(h=1.0, lt=2)
lines(1:nt, st, lwd=2, col="blue")
lines(1:nt, e2, lwd=2, col="red")
legend("topleft", legend=c("mean", "spread", "member"),
       lty=1, col=c("black", "blue", "red"))
