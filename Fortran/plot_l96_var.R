source("confplot_l96_var.R")

title <- paste("Var", "step=", a)
pngname <- "l96_var.png"
w <- 1024
h <-  800

l2 <- readBin("l96_var.dat", "numeric", nc)

png(pngname, w, h)
plot(seq(1, nw * nc, by=nw), l2, type="l", lwd=3,
     main=title, xlab="time", ylab="L2", ylim=c(0, 3),
     cex.lab=2, cex.axis=2, cex.main=3)
abline(h=sqrt(r), lty=2, lwd=2)
abline(h=mean(l2), lty=3, lwd=2)
legend("topleft", cex=2,
       legend=c("analysis",
                paste("mean", format(mean(l2), digits=3)),
                "observation"),
       lty=c(1, 3, 2), lwd=c(3, 2, 2))
dev.off()
