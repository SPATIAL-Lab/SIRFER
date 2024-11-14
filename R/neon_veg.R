source("R/vegFuncs.R")

process.veg("240930_24-251-252.xls")
process.veg("241001_24-252.2.xls")
process.veg("24-252.2-rerun.xls")

veg = read.csv("db/veg.csv")
