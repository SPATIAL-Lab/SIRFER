source("R/vegFuncs.R")


# Five files for processing
process.veg("241118_24-333.1.xls") # noisy spinach d15N
process.veg("241119_24-333.2.xls") # opposite d15N trend 
process.veg("241117_24-251-252.1-RERUN.xls") # crap
process.veg("241114_24-315.2-320.xls") #
process.veg("241115_24-315.1-24-320.2.xls") #
process.veg("241121_24-252.2-RERUN-2.xls") # looks good


# Generate a report
report.veg("manifest_for_D0620240710111247324-24-251.csv")
report.veg("manifest_for_D1120240812104411337-24-315.csv")
report.veg("manifest_for_D1320240807142613008-24-320.csv")


report.veg("manifest_for_D0620240710111247324-24-251.csv")
report.veg("manifest_for_D0620240710111247324-24-251.csv")






