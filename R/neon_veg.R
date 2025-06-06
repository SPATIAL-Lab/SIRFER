source("R/vegFuncs.R")

# Processed
process.veg("241118_24-333.1.xls") # noisy spinach d15N
process.veg("241119_24-333.2.xls") # opposite d15N trend 
process.veg("241117_24-251-252.1-RERUN.xls") # crap
process.veg("241114_24-315.2-320.xls") #
process.veg("241115_24-315.1-24-320.2.xls") #
process.veg("241121_24-252.2-RERUN-2.xls") # looks good
process.veg("241127_24-364.1.xls") # looks good
process.veg("241202_24-364.2.xls") # looks good
process.veg("241204_24-333.1-RERUN.xls")
process.veg("241203_24-379.1.xls")
process.veg("241205_24-379.2.xls")
process.veg("241212_24-379.3-383.xls")
process.veg("250311_24-439.2-240.1.xls")
process.veg("250315_24-439.1.xls")
process.veg("250312_24-440.2.xls")
process.veg("250317_24-440.3-441.1.xls")
## Read this one for low d15N value
View(read.veg("250317_24-440.3-441.1.xls"))
process.veg("250318_24-441.2.xls")
## Read this one for low yield
View(read.veg("250318_24-441.2.xls"))
process.veg("250319_24-441.3.xls")
## Low yield
View(read.veg("250319_24-441.3.xls"))
process.veg("250320_24-441.4-443.xls")
process.veg("250321_24-442.1.xls")
process.veg("250323_24-442.2-441.5.xls")

# Not yet processed

# Reported
report.veg("manifest_for_D1120240812104411337.csv")
report.veg("manifest_for_D0820240611135002668.csv")
report.veg("manifest_for_D0620240710111247324.csv")
report.veg("manifest_for_D0520240913150338770.csv")
report.veg("manifest_for_D0320240801082752226.csv")
report.veg("manifest_for_D1720241010100424502.csv")
report.veg("manifest_for_D0620241014123903736.csv")
report.veg("manifest_for_D1320240807142613008-24-320.csv")
report.veg("manifest_for_D192024120312931971-24-439.csv")
report.veg("manifest_for_D1920241121120601456-24-443.csv")
report.veg("manifest_for_D1320241203120718173-24-441.csv", TRUE)
report.veg("manifest_for_D06202412994102391-24-442.csv")
report.veg("manifest_for_D0120241202120802894-24-440.csv")

# Not yet reported




