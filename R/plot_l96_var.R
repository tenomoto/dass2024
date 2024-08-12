source("config_l96.R")

con <- file("l96_var.dat", "rb")
l2 <- readBin(con, "numeric", nc)
close(con)
title <- paste("4DVar", "step=", a)
plot(seq(1, nt, by=nw), l2, type="l", lwd=2,
     main=title, xlab="time", ylab="l2", ylim=c(0, 3))
abline(h=1, lt=2)
abline(h=mean(l2), lt=3)
legend("topleft", legend=c(
  "analysis", paste("mean", format(mean(l2), digits=3)), "observation"),
  lty=c(1, 3, 2))