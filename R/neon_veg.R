library(readxl)
library(segmented)

# Read file
d = as.data.frame(
  read_xls("data/241001_24-252.2.xls", sheet = "N_conflo.wke"))

# Pick only peak 2
d = d[d$`Peak Nr` == 2,]

# Drop blanks and conditioners
d = d[!(d$`Identifier 1` %in% c("blank tin", "COND")),]

# Add ln(Area)
d$ln_Area = log(d$`Area All N`)

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

# Extract PLRMs
plrm1 = d[d$`Identifier 1` == "UU-CN-3",]
plrm2 = d[d$`Identifier 1` == "UU-CN-2",]

# Calibration equation
cal.s = (9.3 + 4.6) / (mean(plrm1$d15N_lc) - mean(plrm2$d15N_lc)) 
cal.i = -4.6 - mean(plrm2$d15N_lc) * cal.s

# Calibrate all
d$d15N_cal = d$d15N_lc * cal.s + cal.i

# %N calibration
n.s = mean(plrm2$Amount) * 0.0952 / mean(plrm2$`Area All N`)

# Calibrate all
d$Npct = d$`Area All N` * n.s / d$Amount # In one file this field name is `Area All`

# Re-extract all RMs
plrm1 = d[d$`Identifier 1` == "UU-CN-3",]
plrm2 = d[d$`Identifier 1` == "UU-CN-2",]
slrm = d[d$`Identifier 1` == "SPINACH",]

# QC report
data.frame("ID" = c("UU-CN-3", "UU-CN-2", "SPINACH"),
           "d15N_known" = c(9.3, -4.6, -0.4),
           "d15N_cal" = c(mean(plrm1$d15N_cal), mean(plrm2$d15N_cal),
                          mean(slrm$d15N_cal)),
           "d15N_cal.sd" = c(sd(plrm1$d15N_cal), sd(plrm2$d15N_cal),
                          sd(slrm$d15N_cal)),
           "Npct_known" = c(NA, 9.52, 5.95),
           "Npct_meas" = c(mean(plrm1$Npct), mean(plrm2$Npct), mean(slrm$Npct)) * 100,
           "Npct_meas.sd" = c(sd(plrm1$Npct), sd(plrm2$Npct), sd(slrm$Npct)) * 100)

#### Need screening criteria for peak area relative to SLRM, N yield