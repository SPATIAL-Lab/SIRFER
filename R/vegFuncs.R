process.veg = function(fn){
  
  # Prep data
  d = prepVeg(file.path("data", fn))
  
  # Linearity fit
  lin.seg = lin.fit(d)
  
  # Linearity correction
  d.cor = lin.cor(d, lin.seg)
  
  # Calibration
  d.cal = calibrate(d.cor)
  
  # QC report
  d.qc = QC(d.cal)
  
  # Write results
  write.veg(d.qc, fn)
  
}

prepVeg = function(fn){
  library(readxl)
  
  # Read file
  d.N = as.data.frame(read_xls(fn, sheet = "N_conflo.wke"))
  d.C = as.data.frame(read_xls(fn, sheet = "C_conflo.wke"))
  
  # Pick peaks
  d.N = d.N[d.N$`Peak Nr` == 2,]
  d.C = d.C[d.C$`Peak Nr` == 3,]
  
  # Check column names
  names(d.N) = gsub("Area All", "AreaAllN", names(d.N))
  names(d.N) = gsub("%", "PctN", names(d.N))
  names(d.N) = gsub("/", "_", names(d.N))
  names(d.N) = gsub("Start", "StartN", names(d.N))
  names(d.N) = gsub("End", "EndN", names(d.N))

  names(d.C) = gsub("Area All", "AreaAllC", names(d.C))
  names(d.C) = gsub("%", "PctC", names(d.C))
  names(d.C) = gsub("/", "_", names(d.C))
  names(d.C) = gsub("Start", "StartC", names(d.C))
  names(d.C) = gsub("Width", "WidthC", names(d.C))
  
  # Merge
  d = cbind(d.N, d.C[, 8:21])
  
  # Remove WS from field names
  names(d) = gsub(" ", "", names(d))

  # Drop blanks and conditioners
  d = d[!(d$Identifier1 %in% c("blank tin", "COND")),]
  
  # Add ln(Area)
  d$ln_AreaN = log(d$AreaAllN)
  
  return(d)
}

lin.fit = function(d){
  library(segmented)
  
  # Parse out slrms
  slrm = d[d$Identifier1 == "SPINACH",]
  
  # Plot linearity
  opar = par("mar")
  par(mar = c(5, 5, 1, 1))
  plot(slrm$ln_AreaN, slrm$d15N_14N,
       xlab = "ln(AreaN)", ylab = expression("SPINACH  "*delta^15*"N"))
  par(mar = opar)
  
  # Segmented regression
  lin = lm(d15N_14N ~ ln_AreaN, data = slrm)
  lin.seg = segmented(lin, ~ ln_AreaN, psi = 4, 
                      control = seg.control(alpha = c(0.25, 0.75)))

  # Add fit to plot, report r^2
  plot(lin.seg, add = TRUE)  
  cat("Breakpoint =", summary(lin.seg)$psi[2], "\n")
  cat("R2 =", summary(lin.seg)$r.squared, "\n\n")
  
  return(lin.seg)
}

lin.cor = function(d, lin.seg){
  # Correct all data
  d15N_pred = predict(lin.seg, data.frame("ln_AreaN" = d$ln_AreaN))
  
  d$d15N_lc = d$d15N_14N - d15N_pred

  return(d)
}

calibrate = function(d){
  # Extract PLRMs
  plrm1 = d[d$Identifier1 == "UU-CN-3",]
  plrm2 = d[d$Identifier1 == "UU-CN-2",]
  
  # Calibration equations
  calN.s = (9.3 + 4.6) / (mean(plrm1$d15N_lc) - mean(plrm2$d15N_lc)) 
  calN.i = -4.6 - mean(plrm2$d15N_lc) * calN.s
  calC.s = (-12.35 + 28.18) / (mean(plrm1$d13C_12C) -
                                 mean(plrm2$d13C_12C))
  calC.i = -28.18 - mean(plrm2$d13C_12C) * calC.s

  n.s = mean(plrm2$Amount) * 9.52 / mean(plrm2$AreaAllN)
  c.s = mean(plrm2$Amount) * 40.81 / mean(plrm2$AreaAllC)
  
  # Calibrate all
  d$d15N_cal = d$d15N_lc * calN.s + calN.i
  d$d13C_cal = d$d13C_12C * calC.s + calC.i

  d$Npct = d$AreaAllN * n.s / d$Amount
  d$Cpct = d$AreaAllC * c.s / d$Amount
  
  return(d)
}

QC = function(d){
  # Extract all RMs
  plrm1 = d[d$Identifier1 == "UU-CN-3",]
  plrm2 = d[d$Identifier1 == "UU-CN-2",]
  slrm = d[d$Identifier1 == "SPINACH",]
  d = d[!(d$Identifier1 %in% c("UU-CN-3", "UU-CN-2", "SPINACH")),]

  d15N_known = c(9.3, -4.6, -0.4)  
  Npct_known = c(NA, 9.52, 5.95)
  d13C_known = c(-12.35, -28.18, -27.41)
  Cpct_known = c(NA, 40.81, 40.53)
  
  d15N_cal = round(c(mean(plrm1$d15N_cal), mean(plrm2$d15N_cal),
               mean(slrm$d15N_cal)), 1)
  d15N_cal.sd = round(c(sd(plrm1$d15N_cal), sd(plrm2$d15N_cal),
                        sd(slrm$d15N_cal)), 2)
  d13C_cal = round(c(mean(plrm1$d13C_cal), mean(plrm2$d13C_cal),
                     mean(slrm$d13C_cal)), 1)
  d13C_cal.sd = round(c(sd(plrm1$d13C_cal), sd(plrm2$d13C_cal),
                        sd(slrm$d13C_cal)), 2)
  
  Npct_meas = round(c(mean(plrm1$Npct), mean(plrm2$Npct), 
                      mean(slrm$Npct)), 1)
  Npct_meas.sd = round(c(sd(plrm1$Npct), sd(plrm2$Npct), 
                         sd(slrm$Npct)), 1)
  Cpct_meas = round(c(mean(plrm1$Cpct), mean(plrm2$Cpct), 
                      mean(slrm$Cpct)), 1)
  Cpct_meas.sd = round(c(sd(plrm1$Cpct), sd(plrm2$Cpct), 
                         sd(slrm$Cpct)), 1)
  
  d15N.flag = d15N_sd.flag = Npct.flag = Npct_sd.flag =
    d13C.flag = d13C_sd.flag = Cpct.flag = Cpct_sd.flag = rep("", 3)
  
  # QC criteria
  for(i in 1:3){
    if(abs(d15N_cal[i] - d15N_known[i]) > 0.4){
      d15N.flag[i] = "*"
    }
    
    if(!is.na(Npct_known[i]) & abs(Npct_meas[i] - Npct_known[i]) > 0.6){
      Npct.flag[i] = "*"
    }
    
    if(d15N_cal.sd[i] > 0.4){
      d15N_sd.flag[i] = "*"
    }
    
    if(Npct_meas.sd[i] > 0.6){
      Npct_sd.flag[i] = "*"
    }
    
    if(abs(d13C_cal[i] - d13C_known[i]) > 0.4){
      d13C.flag[i] = "*"
    }
    
    if(!is.na(Cpct_known[i]) & abs(Cpct_meas[i] - Cpct_known[i]) > 0.6){
      Cpct.flag[i] = "*"
    }
    
    if(d13C_cal.sd[i] > 0.4){
      d13C_sd.flag[i] = "*"
    }
    
    if(Cpct_meas.sd[i] > 0.6){
      Cpct_sd.flag[i] = "*"
    }
  }
  
  # QC report
  print(data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
             "d15N_known" = as.character(d15N_known),
             "d15N_cal" = paste0(d15N_known, d15N.flag),
             "d15N_cal.sd" = paste0(d15N_cal.sd, d15N_sd.flag),
             "Npct_known" = as.character(Npct_known),
             "Npct_meas" = paste0(Npct_meas, Npct.flag),
             "Npct_meas.sd" = paste0(Npct_meas.sd, Npct_sd.flag)))

  print(data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
                   "d13C_known" = as.character(d13C_known),
                   "d13C_cal" = paste0(d13C_known, d13C.flag),
                   "d13C_cal.sd" = paste0(d13C_cal.sd, d13C_sd.flag),
                   "Cpct_known" = as.character(Cpct_known),
                   "Cpct_meas" = paste0(Cpct_meas, Cpct.flag),
                   "Cpct_meas.sd" = paste0(Cpct_meas.sd, Cpct_sd.flag)))
  
  # Recombine w/ standards
  d = rbind(d, plrm1, plrm2, slrm)
  
  # Add QC flags back to d
  if(Npct_sd.flag[3] == "*" | Cpct_sd.flag[3] == "*"){
    d$cnPercentQF = rep(1)
  }
  else{
    d$cnPercentQF = rep(0)
  }
  
  if(d15N_sd.flag[3] == "*" | d13C_sd.flag[3] == "*"){
    d$cnIsotopeQF = rep(1)
  }
  else{
    d$cnIsotopeQF = rep(0)
  }
  
  if(Npct.flag[3] == "*" | Cpct.flag[3] == "*"){
    d$percentAccuracyQF = rep(1)
  }
  else{
    d$percentAccuracyQF = rep(0)
  }
  
  if(d15N.flag[3] == "*" | d13C.flag[3] == "*"){
    d$isotopeAccuracyQF = rep(1)
  }
  else{
    d$isotopeAccuracyQF = rep(0)
  }
  
  # Sample QC
  d$yieldQF = rep(0)
  d$yieldQF[d$AreaAllN > max(slrm$AreaAllN) * 1.25 |
      d$AreaAllN < min(slrm$AreaAllN) * 0.75 |
        d$AreaAllC > max(slrm$AreaAllC) * 1.25 |
        d$AreaAllC < min(slrm$AreaAllC) * 0.75] = 1
  d$yieldQF[d$Amount * d$Npct < 0.015] = 2
  
  if(sum(d$yieldQF == 2) > 0){
  cat("Small samples\n")
  cat(paste(d$Line[d$yieldQF == 2], 
            d$Identifier1[d$yieldQF == 2], "\n"))
  }
    
  return(d)
}

write.veg = function(d, fn){
  if(!dir.exists("db")){
    dir.create("db")
  }
  
  d$DataFile = rep(fn)
  
  if(file.exists("db/veg.csv")){
    # Read existing data
    veg = read.csv("db/veg.csv")
    
    # Find and remove existing data for these analyses, if any
    v.ui = paste0(veg$Identifier1, veg$TimeCode)
    d.ui = paste0(d$Identifier1, d$TimeCode)

    sup = veg[v.ui %in% d.ui, ]
    veg = veg[!(v.ui %in% d.ui), ]
    
    #####
    # TODO - identify prev runs of sample and assign analysis number
    # TODO - identify if CO2 has been trapped, purge C data, flag
    #####
    
    # Append new data
    veg = rbind(veg, d)
    
    # Write DB and superseded records
    write.csv(veg, "db/veg.csv", row.names = FALSE)
    if(nrow(sup) > 0){
      if(!dir.exists("db/sup")){
        dir.create("db/sup")
      }
      
      fn = file.path("db", "sup", paste0(Sys.Date(), ".csv"))
      for(i in 1:9){
        if(file.exists(fn)){
          fn = file.path("db", "sup", paste0(Sys.Date(), "_", i, ".csv"))
        }
      }
      write.csv(sup, fn, row.names = FALSE)
      cat(nrow(sup), "samples were superseded in the database.")
    }
  }
  else{
    veg = d
    write.csv(veg, "db/veg.csv", row.names = FALSE)
  }
}

read.veg = function(fn, rm.incl = FALSE){
  if(!inherits(fn, "character")){
    stop("Must be one or more run file names")
  }
  
  d = read.csv("db/veg.csv")
  d = d[d$DataFile %in% fn,]
  
  if(!rm.incl){
    d = d[!(d$Identifier1 %in% c("UU-CN-3", "UU-CN-2", "SPINACH")),]
  }
  
  return(d)
}

#####
# TODO: Include option to report QF'd analyses
#####
report.veg = function(fn){
  if(!inherits(fn, "character")){
    stop("Must be a manifest file name")
  }
  
  if(length(fn) > 1){
    stop("Only one manifest can be reported at a time")
  }
  
  # Read manifest
  man = read.csv(file.path("data", fn))
  veg = read.csv("db/veg.csv")
  
  # Samples for manifest
  veg.sub = merge(man, veg, by.x = "SIRFER.ID", by.y = "Identifier1")
  
  if(!all(man$SIRFER.ID %in% veg.sub$SIRFER.ID)){
    miss = man$SIRFER.ID[!(man$SIRFER.ID %in% veg.sub$SIRFER.ID)]
    message(paste("Samples", miss, "not in database"))
  }
  
  # RMs for the manifest
  runs = unique(veg.sub$DataFile)
  rm.sub = veg[veg$Identifier1 == "SPINACH",]
  rm.sub = rm.sub[rm.sub$DataFile %in% runs,]
  rm.sub$rep = rep(0)
  for(i in runs){
    nr = sum(rm.sub$DataFile == i)
    rm.sub$rep[rm.sub$DataFile == i] = 1:nr
  }
  
  # Sample reporting
  veg.out = data.frame("analysisDate" = format(as.Date(veg.sub$TimeCode), "%Y%m%d"),
                       "sampleID" = veg.sub$sampleID,
                       "sampleCode" = veg.sub$sampleCode,
                       "sampleType" = rep("vegetation"),
                       "internalLabID" = veg.sub$SIRFER.ID,
                       "runID" = gsub(".xls", "", veg.sub$DataFile),
                       "acidTreatment" = rep(""),
                       "co2Trapped" = rep("N"),
                       "nitrogenPercent" = veg.sub$Npct,
                       "carbonPercent" = veg.sub$Cpct,
                       "CNratio" = veg.sub$Cpct / veg.sub$Npct,
                       "d15N" = veg.sub$d15N_cal,
                       "d13C" = veg.sub$d13C_cal,
                       "analyticalRepNumber" = rep(1),
                       "cnPercentQF" = veg.sub$cnPercentQF,
                       "cnIsotopeQF" = veg.sub$cnIsotopeQF,
                       "percentAccuracyQF" = veg.sub$percentAccuracyQF,
                       "isotopeAccuracyQF" = veg.sub$isotopeAccuracyQF,
                       "remarks" = rep(""),
                       "testMethod" = rep("NEON_vegIso_SOP v1.0"),
                       "instrument" = rep("Delta Advantage coupled with EA via Conflo III"),
                       "analyzedBy" = rep("schakraborty"),
                       "reviewedBy" = rep("schakraborty"))

  # QC reporting
  ref.out = data.frame("analysisDate" = format(as.Date(rm.sub$TimeCode), "%Y%m%d"),
                       "qaReferenceID" = rep("Spinach"),
                       "internalLabID" = rep(NA),
                       "runID" = gsub(".xls", "", rm.sub$DataFile),
                       "nitrogenPercent" = rm.sub$Npct,
                       "carbonPercent" = rm.sub$Cpct,
                       "CNratio" = rm.sub$Cpct / rm.sub$Npct,
                       "d15N" = rm.sub$d15N_cal,
                       "d13C" = rm.sub$d13C_cal,
                       "analyticalRepNumber" = rm.sub$rep,
                       "cnPercentQF" = rm.sub$cnPercentQF,
                       "cnIsotopeQF" = rm.sub$cnIsotopeQF,
                       "percentAccuracyQF" = rm.sub$percentAccuracyQF,
                       "isotopeAccuracyQF" = rm.sub$isotopeAccuracyQF,
                       "remarks" = rep(""),
                       "testMethod" = rep("NEON_vegIso_SOPv.1.0"),
                       "instrument" = rep("Delta Advantage coupled with EA via Conflo III"),
                       "analyzedBy" = rep("schakraborty"),
                       "reviewedBy" = rep("gjbowen"))

  # Save report files
  dfn = file.path("out", substr(fn, regexec("D", fn)[[1]][1], nchar(fn) - 4))
  qfn = paste0(dfn, "_QA.csv")
  dfn = paste0(dfn, ".csv")
  write.csv(veg.out, dfn, row.names = FALSE)
  write.csv(ref.out, qfn, row.names = FALSE)

  # Summary
  cat(nrow(man), "samples in manifest.\n")
  cat(nrow(veg.out), "rows reported.\n")
  qf = apply(veg.out[, 14:17] != 0, 1, any)
  cat(sum(qf), "rows flagged for quality.\n")
}

