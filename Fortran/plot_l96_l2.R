source("config.R")

title <- paste("4DVar", "step=", a)
fname <- "l2_l96_var.png"
w <- 1024
h <-  800

l2 <- readBin("l2.dat", "numeric", nc)

png(fname, w, h)
plot(seq(1, nw * nc, by=nw), l2, type="l", lwd=3,
     main=title, xlab="time", ylab="l2", ylim=c(0, 3),
     cex.lab=2, cex.axis=2, cex.main=3)
abline(h=1, lty=2, lwd=2)
abline(h=mean(l2), lty=3, lwd=2)
legend("topleft", cex=2,
       legend=c("analysis",
                paste("mean", format(mean(l2), digits=3)),
                "observation"),
       lty=c(1, 3, 2), lwd=c(3, 2, 2))
dev.off()
