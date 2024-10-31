source("R/vegFuncs.R")

# Prep data
d = prepVeg("data/241001_24-252.2.xls")

# Linearity fit
lin.seg = lin.fit(d)

# Linearity correction
d.cor = lin.cor(d, lin.seg)

# Calibration
d.cal = calibrate(d.cor)

# QC report
d.qc = report(d.cal)

