source("R/vegFuncs.R")

process.veg("240930_24-251-252.xls")
process.veg("241001_24-252.2.xls")
process.veg("24-252.2-rerun.xls")






d1 = read.veg("241001_24-252.2.xls")
d2 = read.veg("24-252.2-rerun.xls")

d1.ind = match(d1$Identifier1, d2$Identifier1)
d1.ind = d1.ind[!is.na(d1.ind)]

png("~/252.png", width = 9, height = 4, units = "in", res = 600)
layout(matrix(c(1, 2), ncol = 2))
par(mar = c(5, 5, 1, 1))
plot(d1$d13C_cal[d1.ind], d2$d13C_cal, 
     xlab = expression(delta^13*"C run 1"),
     ylab = expression(delta^13*"C run 2"))
abline(0, 1)
mo = round(mean(d1$d13C_cal[d1.ind] - d2$d13C_cal), 2)
legend("topleft", paste("Mean offset =", mo), bty = "n")

plot(d1$d15N_cal[d1.ind], d2$d15N_cal, 
     xlab = expression(delta^15*"N run 1"),
     ylab = expression(delta^15*"N run 2"))
abline(0, 1)
mo = round(mean(d1$d15N_cal[d1.ind] - d2$d15N_cal), 2)
legend("topleft", paste("Mean offset =", mo), bty = "n")

dev.off()
