source("R/vegFuncs.R")

# Processing ----

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
process.veg("250507_24-364.1-redo.xls")
process.veg("25-065.1.xls")
process.veg("25-065.2.xls")
process.veg("25-065.3.xls")
process.veg("25-065.4.xls")
process.veg("25-066.1.xls")
process.veg("25-066.2.xls")
process.veg("25-066.3.xls")
process.veg("25-021.xls")
process.veg("25-023.xls")
process.veg("25-024.xls")
process.veg("25-068.xls")
process.veg("25-069.1.xls")
process.veg("25-069.2.xls")
process.veg("25-069.3.xls")
process.veg("25-069.4.xls")
process.veg("25-047.1.xls")
process.veg("25-047.2.xls")
process.veg("25-095.1.xls")
process.veg("25-095.2.xls")
process.veg("25-095.3.xls")
process.veg("25-095.4-98.1.xls")
process.veg("25-098.2.xls")
process.veg("25-098.3.xls")
process.veg("25-098.4.xls")
process.veg("25-064.1.xls")
process.veg("25-065.1.xls")
process.veg("25-065.2.xls")
process.veg("25-065.3.xls")
process.veg("25-065.4.xls")
process.veg("25-096.1.xls")
process.veg("25-096.2.xls")
process.veg("25-096.3.xls")
process.veg("25-096.4.xls")
process.veg("25-101.xls")

# Not yet processed

# Reporting ----

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
report.veg("manifest_for_D0520240913150338770-24-364.csv")
report.veg("manifest_for_D032025021291856385-25-065.csv")
report.veg("manifest_for_D172025020392305528-25-066.csv")
report.veg("manifest_for_D1320250106113438993-25-021.csv")
report.veg("manifest_for_D1120250114083929854-25-023.csv")
report.veg("manifest_for_D0120250103140101600-25-024.csv")
report.veg("manifest_for_D012025021193042322-25-069.csv")
report.veg("manifest_for_D1720250128112824458-25-047.csv")
report.veg("manifest_for_D0820250225130525311-25-064.csv")
report.veg("manifest_for_D1720250403134927140-25-095.csv")
report.veg("manifest_for_D1920250207144604785-25-098.csv")
report.veg("manifest_for_D0520250218125157604-25-101.csv")

# Not yet reported
report.veg("manifest_for_D0320250212131213350-25-068.csv") # 10 reanalyses
report.veg("manifest_for_D082025032692651637-25-096.csv") # one reanalysis
