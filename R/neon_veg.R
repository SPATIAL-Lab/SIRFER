source("R/vegFuncs.R")

# Five files for processing
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



# Reported
report.veg("manifest_for_D1120240812104411337.csv")
report.veg("manifest_for_D0820240611135002668.csv")
report.veg("manifest_for_D0620240710111247324.csv")
report.veg("manifest_for_D0520240913150338770.csv")
report.veg("manifest_for_D0320240801082752226.csv")
report.veg("manifest_for_D1720241010100424502.csv")
report.veg("manifest_for_D0620241014123903736.csv")

# Not yet reported

