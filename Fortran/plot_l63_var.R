source("confplot_l63_var.R")

w <- 1024
h <-  800

fname <- paste0("l63_var.dat")
con <- file(fname, "rb")
xt <- matrix(readBin(con, "numeric", ns*nt), nrow=ns)
yo <- matrix(readBin(con, "numeric", ns*nobs), nrow=ns)
xb <- matrix(readBin(con, "numeric", ns*nt), nrow=ns)
close(con)

pngname <- paste0("l63_var.png")
png(pngname, w, h)
tab10.new <- c('#5778a4', '#e49444', '#d1615d', '#85b6b2', '#6a9f58',
               '#e7ca60', '#a87c9f', '#f1a2a9', '#967662', '#b8b0ac')
title <- paste("L63 Var", " ni=", ni, " L2=", format(l2, digits=3))
off <- c(0, 40, 60)
t <- (1:nt) * dt
to <- ((1:nobs) * obs.int) * dt
xt <- xt + off
xb <- xb + off
yo <- yo + off
plot(1, type = "n", xlab="time", ylab="state", main=title,
     xlim=c(min(t), max(t)), ylim=c(-20, 150),
     cex.lab=2, cex.axis=2, cex.main=3)
for (i in 1:ns) {
  lines(t, xt[i,], lwd=3, lty=2, col=tab10.new[i])
  lines(t, xb[i,], lwd=3, col=tab10.new[i])
  points(to, yo[i,], cex=2, pch=4, col=tab10.new[i])
  abline(h=off[i], lwd=2, lty=3, col=tab10.new[i])
}
legend("topright", cex=2, legend=c("mean", "true", "obs"), 
       lty=c(1, 2, NA), pch=c(NA, NA, 4), ncol=3)
dev.off()
