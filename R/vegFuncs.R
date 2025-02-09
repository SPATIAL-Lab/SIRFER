process.veg = function(fn){
  source("R/helpers.R")
  
  # Prep data
  d = prepVeg(file.path("data", fn))
  
  # Correction fit
  corr = corr.fit(d)
  
  # Apply corrections
  d.corr = corr.apply(d, corr)
  
  # Calibration
  d.cal = calibrate(d.corr)
  
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
  
  # Pull comments
  com = unique(d.N$Comment)
  
  # Pick peaks
  d.N = d.N[d.N$`Is Ref _` == 0,]
  d.C = d.C[d.C$`Is Ref _` == 0,]
  
  # Check column names
  names(d.N) = gsub("^Area All$", "AreaAllN", names(d.N))
  names(d.N) = gsub("%", "PctN", names(d.N))
  names(d.N) = gsub("/", "_", names(d.N))
  names(d.N) = gsub("^Start$", "StartN", names(d.N))
  names(d.N) = gsub("^End$", "EndN", names(d.N))

  names(d.C) = gsub("^Area All$", "AreaAllC", names(d.C))
  names(d.C) = gsub("%", "PctC", names(d.C))
  names(d.C) = gsub("/", "_", names(d.C))
  names(d.C) = gsub("^Start$", "StartC", names(d.C))
  names(d.C) = gsub("^Width$", "WidthC", names(d.C))
  
  # Merge
  d = merge(d.N, d.C[, c(1, 8:21)], by = "Line", all = TRUE)
  
  # Remove WS from field names
  names(d) = gsub(" ", "", names(d))
  
  # Catch missing IDs
  for(i in seq_along(d$Line)){
    if(is.na(d$Identifier1[i])){
      d[i, 2:5] = d.C[match(d$Line[i], d.C$Line), 2:5]
    }
  }
  
  # Blank warning
  if("blank tin" %in% d$Identifier1){
    cat("Blank peak detected\n\n")
  }

  # Drop blanks and conditioners
  d = d[!(d$Identifier1 %in% c("blank tin", "COND")),]
  
  # Identify CO2 trapping
  jn = trimws(com[grep("Job", com)])
  if(substr(jn, nchar(jn) - 3, nchar(jn)) %in% c("trap", "Trap", "TRAP")){
    d$Trapping = rep("Y")
    d[, 20:33] = rep(NA)
    cat("CO2 trapping detected\n\n")
  } else{
    d$Trapping = rep("N")
  }

  # Identify missing peaks
  d.missing = d[0, ]
  for(i in seq_along(d$Line)){
    start.m = mean(d$StartN, na.rm = TRUE)
    start.sd = sd(d$StartN, na.rm = TRUE)
    if(is.na(d$StartN[i]) | abs(d$StartN[i] - start.m) > 5 * start.sd){
      d.missing = rbind(d.missing, d[i, ])
    }
  }
  d = d[!(d$Line %in% d.missing$Line), ]

  if(d$Trapping[1] == "N"){
    for(i in seq_along(d$Line)){
      start.m = mean(d$StartC, na.rm = TRUE)
      start.sd = sd(d$StartC, na.rm = TRUE)
      if(is.na(d$StartC[i]) | abs(d$StartC[i] - start.m) > 5 * start.sd){
        d.missing = rbind(d.missing, d[i, ])
      }
    }
  }
  d = d[!(d$Line %in% d.missing$Line), ]
  
  if(nrow(d.missing) > 0){
    cat("One or more missing peaks\n")
    cat(paste(d.missing$Line,
              d.missing$Identifier1, "\n", sep = "\t"), "\n")
  }
  
  # Add ln(Area)
  d$ln_AreaN = log(d$AreaAllN)
  d$ln_AreaC = log(d$AreaAllC)
  
  return(d)
}

corr.fit = function(d){
  opar = par("mar")
  on.exit(par(mar = opar))
  par(mar = c(5, 5, 3, 1))

  # Parse out RMs
  plrm1 = d[d$Identifier1 == "UU-CN-3",]
  plrm2 = d[d$Identifier1 == "UU-CN-2",]
  slrm = d[d$Identifier1 == "SPINACH",]

  # Check for drift using plrms
  ## Nitrogen
  drift.n1 = lm(d15N_14N ~ Line, plrm1)
  drift.n2 = lm(d15N_14N ~ Line, plrm2)
  slope.n1 = summary(drift.n1)$coeff[2, 1]
  slope.n2 = summary(drift.n2)$coeff[2, 1]
  if(summary(drift.n1)$adj.r.squared > 0.3 &
     summary(drift.n2)$adj.r.squared > 0.3 &
     identical(sign(slope.n1), sign(slope.n2))){
    slope.n = mean(c(slope.n1, slope.n2))
    cat("Nitrogen drift slope of", round(slope.n, 4), "applied.\n")
    cat("PLRM1 SD before:", round(sd(plrm1$d15N_14N), 2), "and after:",
        round(sd(plrm1$d15N_14N - plrm1$Line * slope.n), 2), "\n")
    cat("PLRM2 SD before:", round(sd(plrm2$d15N_14N), 2), "and after:",
        round(sd(plrm2$d15N_14N - plrm2$Line * slope.n), 2), "\n\n")
  } else{
    slope.n = 0
    cat("No nitrogen drift applied\n\n")
  }
  
  plot(plrm1$Line, plrm1$d15N_14N - mean(plrm1$d15N_14N), xlim = range(d$Line),
       ylim = range(c(plrm1$d15N_14N - mean(plrm1$d15N_14N), 
                      plrm2$d15N_14N - mean(plrm2$d15N_14N))),
       main = "N Drift", xlab = "Line", ylab = expression(Delta*delta^{15}*"N"),
       pch = 21, bg = 2)
  points(plrm2$Line, plrm2$d15N_14N - mean(plrm2$d15N_14N), pch = 21, bg = 3)
  abline(0, slope.n)
  
  ## Carbon
  if(d$Trapping[1] == "N"){
    drift.c1 = lm(d13C_12C ~ Line, plrm1)
    drift.c2 = lm(d13C_12C ~ Line, plrm2)
    slope.c1 = summary(drift.c1)$coeff[2, 1]
    slope.c2 = summary(drift.c2)$coeff[2, 1]
    if(summary(drift.c1)$adj.r.squared > 0.3 &
       summary(drift.c2)$adj.r.squared > 0.3 &
       identical(sign(slope.c1), sign(slope.c2))){
      slope.c = mean(c(slope.c1, slope.c2))
      cat("Carbon drift slope of", round(slope.c, 4), "applied.\n")
      cat("PLRM1 SD before:", round(sd(plrm1$d13C_12C), 2), "and after:",
          round(sd(plrm1$d13C_12C - plrm1$Line * slope.c), 2), "\n")
      cat("PLRM2 SD before:", round(sd(plrm2$d13C_12C), 2), "and after:",
          round(sd(plrm2$d13C_12C - plrm2$Line * slope.c), 2), "\n\n")
    } else{
      slope.c = 0
      cat("No carbon drift correction applied.\n\n")
    }
    
    plot(plrm1$Line, plrm1$d13C_12C - mean(plrm1$d13C_12C), xlim = range(d$Line),
         ylim = range(c(plrm1$d13C_12C - mean(plrm1$d13C_12C), 
                        plrm2$d13C_12C - mean(plrm2$d13C_12C))),
         main = "C Drift", xlab = "Line", ylab = expression(Delta*delta^{13}*"C"),
         pch = 21, bg = 2)
    points(plrm2$Line, plrm2$d13C_12C - mean(plrm2$d13C_12C), pch = 21, bg = 3)
    abline(0, slope.n)
  } else{
    slope.c = 0
  }

  # Apply C and N drift correction to slrms
  slrm$d15N_14N = slrm$d15N_14N - slrm$Line * slope.n
  slrm$d13C_12C = slrm$d13C_12C - slrm$Line * slope.c
  
  # Plot linearity N
  plot(slrm$ln_AreaN, slrm$d15N_14N, main = "N linearity", pch = 21, bg = 2,
       xlab = "ln(AreaN)", ylab = expression("SPINACH  "*delta^15*"N"))

  # Segmented regression
  lin.n = fit.piece(slrm$ln_AreaN, slrm$d15N_14N)
  
  # Add fit to plot, report r^2
  if(lin.n$m != 0){
    lines.piece(lin.n)
    cat("N Breakpoint =", round(lin.n$psi, 2), "\n")
    cat("N R2 =", round(lin.n$r2, 2), "\n\n")
  } else{
    cat("No nitrogen linearity correction applied.\n\n")
  }
  
  if(d$Trapping[1] == "N"){
    # Plot linearity C
    plot(slrm$ln_AreaC, slrm$d13C_12C, main = "C linearity", pch = 21, bg = 2,
         xlab = "ln(AreaC)", ylab = expression("SPINACH  "*delta^13*"C"))
    
    # Segmented regression
    lin.c = fit.piece(slrm$ln_AreaC, slrm$d13C_12C)
    
    # Add fit to plot, report r^2
    if(lin.c$m != 0){
      lines.piece(lin.c)
      cat("C Breakpoint =", lin.c$psi, "\n")
      cat("C R2 =", lin.c$r2, "\n\n")
    } else{
      cat("No carbon linearity correction applied.\n\n")
    }
  } else{
    lin.c = piece(0, 0, 0, 0, 0)
  }
  
  return(list("slope.n" = slope.n, "slope.c" = slope.c, 
              "lin.n" = lin.n, "lin.c" = lin.c))
}

corr.apply = function(d, corr){
  # Apply drift corrections
  d$d15N_dc = d$d15N_14N - d$Line * corr$slope.n
  d$d13C_dc = d$d13C_12C - d$Line * corr$slope.c
  
  # Apply linearity correction
  d15N_pred = predict.piece(corr$lin.n, d$ln_AreaN)
  d$d15N_lc = d$d15N_dc - d15N_pred

  d13C_pred = predict.piece(corr$lin.c, d$ln_AreaC)
  d$d13C_lc = d$d13C_dc - d13C_pred
  
  return(d)
}

calibrate = function(d){
  # Extract PLRMs
  plrm1 = d[d$Identifier1 == "UU-CN-3",]
  plrm2 = d[d$Identifier1 == "UU-CN-2",]
  
  # Calibration equations
  calN.s = (9.3 + 4.6) / (mean(plrm1$d15N_lc) - mean(plrm2$d15N_lc)) 
  calN.i = -4.6 - mean(plrm2$d15N_lc) * calN.s
  calC.s = (-12.35 + 28.18) / (mean(plrm1$d13C_lc) -
                                 mean(plrm2$d13C_lc))
  calC.i = -28.18 - mean(plrm2$d13C_lc) * calC.s

  n.s = mean(plrm2$Amount) * 9.52 / mean(plrm2$AreaAllN)
  c.s = mean(plrm2$Amount) * 40.81 / mean(plrm2$AreaAllC)
  
  # Calibrate all
  d$d15N_cal = d$d15N_lc * calN.s + calN.i
  d$d13C_cal = d$d13C_lc * calC.s + calC.i

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
  
  opar = par("mar")
  on.exit(par(mar = opar))
  par(mar = c(5, 5, 1, 1))
  if(d$Trapping[1] == "N"){
    plot(d$d13C_cal, d$d15N_cal, xlab = expression("Sample "*delta^13*"C"),
         ylab = expression("Sample "*delta^15*"N"), pch = 21, bg = "grey60")
  } else{
    plot(density(d$d15N_cal, na.rm = TRUE), 
         xlab = expression("Sample "*delta^15*"N"), ylab = "",
         main = "", lwd = 2)
  }

  d15N_known = c(9.3, -4.6, -0.4)  
  Npct_known = c(NA, 9.52, 5.95)
  d13C_known = c(-12.35, -28.18, -27.41)
  Cpct_known = c(NA, 40.81, 40.53)
  
  d15N_cal = c(mean(plrm1$d15N_cal), mean(plrm2$d15N_cal),
               mean(slrm$d15N_cal))
  d15N_cal.sd = c(sd(plrm1$d15N_cal), sd(plrm2$d15N_cal),
                        sd(slrm$d15N_cal))
  Npct_meas = c(mean(plrm1$Npct), mean(plrm2$Npct), 
                      mean(slrm$Npct))
  Npct_meas.sd = c(sd(plrm1$Npct), sd(plrm2$Npct), 
                         sd(slrm$Npct))
  
  d13C_cal = c(mean(plrm1$d13C_cal), mean(plrm2$d13C_cal),
                     mean(slrm$d13C_cal))
  d13C_cal.sd = c(sd(plrm1$d13C_cal), sd(plrm2$d13C_cal),
                        sd(slrm$d13C_cal))
  Cpct_meas = c(mean(plrm1$Cpct), mean(plrm2$Cpct), 
                      mean(slrm$Cpct))
  Cpct_meas.sd = c(sd(plrm1$Cpct), sd(plrm2$Cpct), 
                         sd(slrm$Cpct))
  
  d15N.flag = d15N_sd.flag = Npct.flag = Npct_sd.flag =
    d13C.flag = d13C_sd.flag = Cpct.flag = Cpct_sd.flag = rep("", 3)
  
  # QC criteria
  if(abs(d15N_cal[3] - d15N_known[3]) > 0.4){
    d15N.flag[3] = "*"
  }
  
  if(abs(Npct_meas[3] - Npct_known[3]) > 0.6){
    Npct.flag[3] = "*"
  }
  
  if(d15N_cal.sd[3] > 0.4){
    d15N_sd.flag[3] = "*"
  }
  
  if(Npct_meas.sd[3] > 0.6){
    Npct_sd.flag[3] = "*"
  }
  
  if(d$Trapping[1] == "N"){
    if(abs(d13C_cal[3] - d13C_known[3]) > 0.4){
      d13C.flag[3] = "*"
    }
    
    if(abs(Cpct_meas[3] - Cpct_known[3]) > 0.6){
      Cpct.flag[3] = "*"
    }
    
    if(d13C_cal.sd[3] > 0.4){
      d13C_sd.flag[3] = "*"
    }
    
    if(Cpct_meas.sd[3] > 0.6){
      Cpct_sd.flag[3] = "*"
    }
  }
  
  # QC report
  print(data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
             "d15N_known" = as.character(d15N_known),
             "d15N_cal" = paste0(round(d15N_cal, 2), d15N.flag),
             "d15N_cal.sd" = paste0(round(d15N_cal.sd, 2), d15N_sd.flag),
             "Npct_known" = as.character(Npct_known),
             "Npct_meas" = paste0(round(Npct_meas, 2), Npct.flag),
             "Npct_meas.sd" = paste0(round(Npct_meas.sd, 2), Npct_sd.flag)))

  print(data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
                   "d13C_known" = as.character(d13C_known),
                   "d13C_cal" = paste0(round(d13C_cal, 2), d13C.flag),
                   "d13C_cal.sd" = paste0(round(d13C_cal.sd, 2), d13C_sd.flag),
                   "Cpct_known" = as.character(Cpct_known),
                   "Cpct_meas" = paste0(round(Cpct_meas, 2), Cpct.flag),
                   "Cpct_meas.sd" = paste0(round(Cpct_meas.sd, 2), Cpct_sd.flag)))
  
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
  if(d$Trapping[1] == "N"){
    d$yieldQF[d$AreaAllN > max(slrm$AreaAllN) * 1.25 |
                d$AreaAllN < min(slrm$AreaAllN) * 0.75 |
                d$AreaAllC > max(slrm$AreaAllC) * 1.25 |
                d$AreaAllC < min(slrm$AreaAllC) * 0.75] = 1
  } else{
    d$yieldQF[d$AreaAllN > max(slrm$AreaAllN) * 1.25 |
                d$AreaAllN < min(slrm$AreaAllN) * 0.75] = 1
  }
  
  if(sum(d$yieldQF == 1) > 0){
    cat("C or N yield out of range\n")
    cat(paste(d$Line[d$yieldQF == 1], 
              d$Identifier1[d$yieldQF == 1], "\n", sep = "\t"), "\n")
  }
  
  d$yieldQF[d$Amount * d$Npct < 0.015] = 2
  
  if(sum(d$yieldQF == 2) > 0){
  cat("N yield below limit\n")
  cat(paste(d$Line[d$yieldQF == 2], 
            d$Identifier1[d$yieldQF == 2], "\n", sep = "\t"), "\n")
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

report.veg = function(fn, flagged = FALSE){
  opar = par("mar")
  on.exit(par(mar = opar))
  par(mar = c(5, 5, 3, 1))
  
  if(!inherits(fn, "character")){
    stop("fn must be a manifest file name")
  }
  
  if(length(fn) > 1){
    stop("Only one manifest can be reported at a time")
  }
  
  if(!inherits(flagged, "logical")){
    stop("flagged must be TRUE/FALSE")
  }
  
  # Read manifest
  man = read.csv(file.path("data", fn))
  veg = read.csv("db/veg.csv")
  
  # Protect against inconsistent manifests
  man = man[c("SIRFER.ID", "domainID", "sampleID", "sampleCode")]
  if(ncol(man) != 4){
    stop("Missing columns in manifest")
  }
  
  # Samples for manifest
  veg.sub = merge(man, veg, by.x = "SIRFER.ID", by.y = "Identifier1")
  
  if(!all(man$SIRFER.ID %in% veg.sub$SIRFER.ID)){
    miss = man$SIRFER.ID[!(man$SIRFER.ID %in% veg.sub$SIRFER.ID)]
    message(paste("Sample", miss, "not in database\n"))
  }
  
  # Remove QC flagged if requested
  # Otherwise report flagged samples w/o unflagged analysis
  qf = apply(veg.sub[c("cnPercentQF", "cnIsotopeQF", "percentAccuracyQF", 
                       "isotopeAccuracyQF", "yieldQF")] == 1, 1, any)
  il = nrow(veg.sub)
  if(!flagged){
    veg.sub = veg.sub[!qf, ]
    cat(il - nrow(veg.sub), "rows removed based on QC flags\n")
  } else{
    veg.temp = veg.sub[!qf, ]
    veg.add = veg.sub[!(veg.sub$sampleCode %in% veg.temp$sampleCode),]
    veg.sub = rbind(veg.temp, veg.add)
    cat(il - nrow(veg.sub), "rows removed based on QC flags\n")
  }
  
  # Remove low N data if a CO2 trapped value is available
  qf = veg.sub["yieldQF"] == 2
  if(!flagged){
    veg.sub$d15N_cal[qf] = veg.sub$Npct[qf] = NA
    cat(sum(is.na(veg.sub$d15N_cal)), "N values removed for low yield\n")
  } else{
    veg.lowN = veg.sub$sampleCode[qf]
    veg.trap = veg.sub$sampleCode[veg.sub$Trapping == "Y"]
    veg.rm = veg.lowN[veg.lowN %in% veg.trap]
    veg.sub$d15N_cal[qf & veg.sub$sampleCode %in% veg.rm] = 
      veg.sub$Npct[qf & veg.sub$sampleCode %in% veg.rm] = NA
    cat(sum(is.na(veg.sub$d15N_cal)), "N values removed for low yield\n")
  }
  
  # Catch zero
  if(nrow(veg.sub) == 0){
    stop("\nNo samples in report")
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
                       "co2Trapped" = veg.sub$Trapping,
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
                       "testMethod" = rep("NEON_vegIso_SOPv1.0"),
                       "instrument" = rep("Carlo Erba 1110 Elemental Analyzer with Costech Zero Blank Autosampler coupled to Thermo Delta Plus Advantage IRMS with Conflo III Interface"),
                       "analyzedBy" = rep("schakraborty"),
                       "reviewedBy" = rep("gjbowen"))

  plot(veg.out$d13C, veg.out$d15N, main = "Sample isotopes", 
       xlab = expression(delta^{13}*"C"), ylab = expression(delta^{15}*"N"),
       pch = 21, bg = 2)
  
  plot(veg.out$carbonPercent, veg.out$nitrogenPercent, 
       main = "Sample composition", xlab = "wt% C", ylab = "wt% N",
       pch = 21, bg = 2)
  
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
                       "testMethod" = rep("NEON_vegIso_SOP.v1.0"),
                       "instrument" = rep("Carlo Erba 1110 Elemental Analyzer with Costech Zero Blank Autosampler coupled to Thermo Delta Plus Advantage IRMS with Conflo III Interface"),
                       "analyzedBy" = rep("schakraborty"),
                       "reviewedBy" = rep("gjbowen"))

  ## Plot SPINACH isotopes 
  x1 = -27.41 + c(-1, 1, 1, -1) * 0.2
  y1 = -0.4 + c(-1, -1, 1, 1) * 0.2
  x2 = -27.41 + c(-1, 1, 1, -1) * 0.4
  y2 = -0.4 + c(-1, -1, 1, 1) * 0.4
  plot(ref.out$d13C, ref.out$d15N, main = "SPINACH isotopes", 
       xlab = expression(delta^{13}*"C"), ylab = expression(delta^{15}*"N"), 
       xlim = range(c(x2, ref.out$d13C)), 
       ylim = range(c(y2, ref.out$d15N)),
       pch = 21, bg = 2)
  polygon(x1, y1, border = 2, lwd = 2)
  polygon(x2, y2, border = 2, lwd = 2)
  
  ## Plot SPINACH wt% 
  x1 = 40.53 + c(-1, 1, 1, -1) * 0.3
  y1 = 5.95 + c(-1, -1, 1, 1) * 0.3
  x2 = 40.53 + c(-1, 1, 1, -1) * 0.6
  y2 = 5.95 + c(-1, -1, 1, 1) * 0.6
  plot(ref.out$carbonPercent, ref.out$nitrogenPercent, 
       main = "SPINACH composition", xlab = "wt% C", ylab = "wt% N",
       xlim = range(c(x2, ref.out$carbonPercent)),
       ylim = range(c(y2, ref.out$nitrogenPercent)),
       pch = 21, bg = 2)
  polygon(x1, y1, border = 2, lwd = 2)
  polygon(x2, y2, border = 2, lwd = 2)
  
  # Save report files
  jobnum = substr(veg.out$internalLabID[1], 1, 6)
  fname = paste0(substr(fn, regexec("D", fn)[[1]][1], nchar(fn) - 4), 
                 "-", jobnum)
  dfn = file.path("out", fname)
  qfn = paste0(dfn, "_QA.csv")
  dfn = paste0(dfn, ".csv")
  write.csv(veg.out, dfn, row.names = FALSE)
  write.csv(ref.out, qfn, row.names = FALSE)

  # Summary
  cat(nrow(man), "samples in manifest.\n")
  cat(length(unique(veg.out$sampleID)), "samples reported.\n")
}

