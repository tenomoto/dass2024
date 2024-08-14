source("confplot_l96_ens.R")

w <- 1024
h <-  800

fname <- paste0("l96_", fil, ".dat")
e <- matrix(readBin(fname, "numeric", nt*3), nrow=nt)

pngname <- paste0("l96_", fil, ".png")
png(pngname, w, h)
title <- paste(fil, "ne=", ne, "infl=", infl, "loc=", c.loc)
plot(1:nt, e[,1], lwd=3,
     main=title, xlab="time", ylab="L2", type="l", ylim=c(0, 3),
     cex.lab=2, cex.axis=2, cex.main=3
)
abline(h=sqrt(r), lt=2)
lines(1:nt, e[,3], lwd=3, col="blue")
lines(1:nt, e[,2], lwd=3, col="red")
legend("topleft", cex=2, legend=c("mean", "stddev", "member"),
       lty=1, lwd=3, col=c("black", "blue", "red"))
dev.off()
