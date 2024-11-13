prepVeg = function(fn){
  library(readxl)
  
  # Read file
  d.N = as.data.frame(read_xls(fn, sheet = "N_conflo.wke"))
  d.C = as.data.frame(read_xls(fn, sheet = "C_conflo.wke"))
  
  # Pick peaks
  d.N = d.N[d.N$`Peak Nr` == 2,]
  d.C = d.C[d.C$`Peak Nr` == 3,]
  
  # Merge
  d = cbind(d.N, d.C[, 8:21])
  
  # Drop blanks and conditioners
  d = d[!(d$`Identifier 1` %in% c("blank tin", "COND")),]
  
  # Add ln(Area)
  d$ln_AreaN = log(d$`Area All N`)
  
  return(d)
}

lin.fit = function(d){
  library(segmented)
  
  # Parse out slrms
  slrm = d[d$`Identifier 1` == "SPINACH",]
  
  # Plot linearity
  plot(slrm$ln_AreaN, slrm$`d 15N/14N`)
  
  # Segmented regression
  lin = lm(`d 15N/14N` ~ ln_AreaN, data = slrm)
  lin.seg = segmented(lin, seg.Z = ~ ln_AreaN, psi = 4)
  
  # Add fit to plot, report r^2
  abline(lin.seg)
  summary(lin.seg)$r.squared
  
  return(lin.seg)
}

lin.cor = function(d, lin.seg){
  # Correct all data
  d$d15N_lc = rep(0)
  for(i in seq_along(d$Line)){
    if(d$ln_AreaN[i] < lin.seg$psi[2]){
      d$d15N_lc[i] = d$`d 15N/14N`[i] + 
        (lin.seg$psi[2] - d$ln_AreaN[i]) * lin.seg$coefficients[2]
    }
    else{
      d$d15N_lc[i] = d$`d 15N/14N`[i]
    }
  }
  
  return(d)
}

calibrate = function(d){
  # Extract PLRMs
  plrm1 = d[d$`Identifier 1` == "UU-CN-3",]
  plrm2 = d[d$`Identifier 1` == "UU-CN-2",]
  
  # Calibration equations
  calN.s = (9.3 + 4.6) / (mean(plrm1$d15N_lc) - mean(plrm2$d15N_lc)) 
  calN.i = -4.6 - mean(plrm2$d15N_lc) * calN.s
  calC.s = (-12.35 + 28.18) / (mean(plrm1$`d 13C/12C`) -
                                 mean(plrm2$`d 13C/12C`))
  calC.i = -28.18 - mean(plrm2$`d 13C/12C`) * calC.s
  
  # Calibrate all
  d$d15N_cal = d$d15N_lc * calN.s + calN.i
  d$d13C_cal = d$`d 13C/12C` * calC.s + calC.i
  
  # Find N area column
#  n.names = c("Area All N", "Area All")
#  n.c = match(n.names[n.names %in% names(d)], names(d))
  
  # %N calibration
#  n.s = mean(plrm2$Amount) * 0.0952 / mean(plrm2[, n.c])
  n.s = mean(plrm2$Amount) * 9.52 / mean(plrm2$`Area All N`)
  c.s = mean(plrm2$Amount) * 40.81 / mean(plrm2$`Area All`)
  
  # Calibrate all
  d$Npct = d$`Area All N` * n.s / d$Amount
  d$Cpct = d$`Area All` * c.s / d$Amount
  
  return(d)
}

QC = function(d){
  # Extract all RMs
  plrm1 = d[d$`Identifier 1` == "UU-CN-3",]
  plrm2 = d[d$`Identifier 1` == "UU-CN-2",]
  slrm = d[d$`Identifier 1` == "SPINACH",]
  d = d[!(d$`Identifier 1` %in% c("UU-CN-3", "UU-CN-2", "SPINACH")),]

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
  
  # Add QC flags back to d
  if(Npct_sd.flag[3] == "*" | Cpct_sd.flag[3] == "*"){
    d$cnPercentQF = rep(1)
  } else{
    d$cnPercentQF = rep(0)
  }
  
  if(d15N_sd.flag[3] == "*" | d13C_sd.flag[3] == "*"){
    d$cnIsotopeQF = rep(1)
  } else{
    d$cnIsotopeQF = rep(0)
  }
  
  if(Npct.flag[3] == "*" | Cpct.flag[3] == "*"){
    d$percentAccuracyQF = rep(1)
  } else{
    d$percentAccuracyQF = rep(0)
  }
  
  if(d15N.flag[3] == "*" | d13C.flag[3] == "*"){
    d$isotopeAccuracyQF = rep(1)
  } else{
    d$isotopeAccuracyQF = rep(0)
  }
  
  # Sample QC - need distinct flag for low N yield
  d$yieldQF = rep(0)
  d$yieldQF[d$`Area All N` > max(slrm$`Area All N`) * 1.25 |
      d$`Area All N` < min(slrm$`Area All N`) * 0.75 |
        d$`Area All` > max(slrm$`Area All`) * 1.25 |
        d$`Area All` < min(slrm$`Area All`) * 0.75] = 1
  d$yieldQF[d$Amount * d$Npct < 0.015] = 2
  
  if(sum(d$yieldQF == 2) > 0){
  cat("Small samples\n")
  cat(paste(d$Line[d$yieldQF == 2], 
            d$`Identifier 1`[d$yieldQF == 2], "\n"))
  }

  return(d)
}

# Read manifest - where to get these, how to identify?
mf = list.files("data", "manifest")
man = read.csv("data/manifest_for_D0620240710111247324.csv")
