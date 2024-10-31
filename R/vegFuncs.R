prepVeg = function(fn){
  library(readxl)
  
  # Read file
  d = as.data.frame(
    read_xls(fn, sheet = "N_conflo.wke"))
  
  # Pick only peak 2
  d = d[d$`Peak Nr` == 2,]
  
  # Drop blanks and conditioners
  d = d[!(d$`Identifier 1` %in% c("blank tin", "COND")),]
  
  # Add ln(Area)
  d$ln_Area = log(d$`Area All N`)
  
  return(d)
}

lin.fit = function(d){
  library(segmented)
  
  # Parse out slrms
  slrm = d[d$`Identifier 1` == "SPINACH",]
  
  # Plot linearity
  plot(slrm$ln_Area, slrm$`d 15N/14N`)
  
  # Segmented regression
  lin = lm(`d 15N/14N` ~ ln_Area, data = slrm)
  lin.seg = segmented(lin, seg.Z = ~ ln_Area, psi = 4)
  
  # Add fit to plot, report r^2
  abline(lin.seg)
  summary(lin.seg)$r.squared
  
  return(lin.seg)
}

lin.cor = function(d, lin.seg){
  # Correct all data
  d$d15N_lc = rep(0)
  for(i in seq_along(d$Line)){
    if(d$ln_Area[i] < lin.seg$psi[2]){
      d$d15N_lc[i] = d$`d 15N/14N`[i] + 
        (lin.seg$psi[2] - d$ln_Area[i]) * lin.seg$coefficients[2]
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
  
  # Calibration equation
  cal.s = (9.3 + 4.6) / (mean(plrm1$d15N_lc) - mean(plrm2$d15N_lc)) 
  cal.i = -4.6 - mean(plrm2$d15N_lc) * cal.s
  
  # Calibrate all
  d$d15N_cal = d$d15N_lc * cal.s + cal.i
  
  # Find N area column
  n.names = c("Area All N", "Area All")
  n.c = match(n.names[n.names %in% names(d)], names(d))
  
  # %N calibration
  n.s = mean(plrm2$Amount) * 0.0952 / mean(plrm2[, n.c])
  
  # Calibrate all
  d$Npct = d[, n.c] * n.s / d$Amount
  
  return(d)
}

report = function(d){
  # Extract all RMs
  plrm1 = d[d$`Identifier 1` == "UU-CN-3",]
  plrm2 = d[d$`Identifier 1` == "UU-CN-2",]
  slrm = d[d$`Identifier 1` == "SPINACH",]
  d = d[!(d$`Identifier 1` %in% c("UU-CN-3", "UU-CN-2", "SPINACH")),]

  d15N_known = c(9.3, -4.6, -0.4)  
  d15N_cal = round(c(mean(plrm1$d15N_cal), mean(plrm2$d15N_cal),
               mean(slrm$d15N_cal)), 1)
  d15N_cal.sd = round(c(sd(plrm1$d15N_cal), sd(plrm2$d15N_cal),
                        sd(slrm$d15N_cal)), 2)
  Npct_known = c(NA, 9.52, 5.95)
  Npct_meas = round(c(mean(plrm1$Npct), mean(plrm2$Npct), mean(slrm$Npct)) * 
                      100, 1)
  Npct_meas.sd = round(c(sd(plrm1$Npct), sd(plrm2$Npct), sd(slrm$Npct)) * 
                         100, 1)
  
  d15N.flag = d15N_sd.flag = Npct.flag = Npct_sd.flag = rep("", 3)
  
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
  }
  
  # QC report
  data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
             "d15N_known" = as.character(d15N_known),
             "d15N_cal" = paste0(d15N_known, d15N.flag),
             "d15N_cal.sd" = paste0(d15N_cal.sd, d15N_sd.flag),
             "Npct_known" = as.character(Npct_known),
             "Npct_meas" = paste0(Npct_meas, Npct.flag),
             "Npct_meas.sd" = paste0(Npct_meas.sd, Npct_sd.flag))
  
  # Find N area column
  n.names = c("Area All N", "Area All")
  n.c = match(n.names[n.names %in% names(d)], names(d))
  
  # Sample QC
  d$QF = rep(0)
  d$QF[d[, n.c] > max(slrm[, n.c]) * 1.25 |
      d[, n.c] < min(slrm[, n.c]) * 0.75 |
      d$Amount * d$Npct < 0.015] = 1
  
  return(d)
}
