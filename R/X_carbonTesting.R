cplot = function(fn){
  source("R/vegFuncs.R")
  d = prepVeg(file.path("data", fn))
  
  sirm = d[d$Identifier1 == "SPINACH",]
  pirms = d[d$Identifier1 %in% c("UU-CN-2", "UU-CN-3"),]
  sam = d[!(d$Line %in% c(sirm$Line, pirms$Line)),]
  
  plot(log(sirm$AreaAllC), sirm$d13C_12C, pch = 21, bg = "grey", cex = 2)
  points(log(sam$AreaAllC), rep(par("usr")[3], nrow(sam)), pch = 21, 
         bg = "blue", cex = 2, xpd = TRUE)
  points(log(pirms$AreaAllC), rep(par("usr")[3], nrow(pirms)), pch = 21, 
         bg = "red", cex = 2, xpd = TRUE)
}

cplot("241127_24-364.1.xls")
cplot("241114_24-315.2-320.xls")
