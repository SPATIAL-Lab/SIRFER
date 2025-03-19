library(openxlsx)

# All results files
lf = list.files("out/", ".xlsx", full.names = TRUE)

# Find the ones with our samples
rf = character()
for(i in lf){
  d = read.xlsx(i)
  j = grep("FE01*", d$Identifier_1)
  if(length(j) != 0){
    rf = append(rf, i)
  }
}

# Read and compile data
d = read.xlsx(rf[1], sheet = "All")
d$run = rep(rf[1])
for(i in 2:length(rf)){
  d2 = read.xlsx(rf[i], sheet = "All")
  d2$run = rep(rf[i])
  d = rbind(d, d2)
}

# Find and fix a few typos
sort(unique(d$Identifier_1))
d$Identifier_1 = gsub("B00-11-CT", "BOO-11-CT", d$Identifier_1)
d$Identifier_1 = gsub("BOO-II*", "BOO-11", d$Identifier_1)
d$Identifier_1 = gsub("CARRARA", "CM", d$Identifier_1)
d$Identifier_1 = gsub("Carrara", "CM", d$Identifier_1)
d$Identifier_1 = gsub("MAR", "Marble", d$Identifier_1)
d = d[d$Identifier_1 != "LSVEC",]

# Split the data
stds = c("CM", "CO-8", "Marble")
rm = d[d$Identifier_1 %in% stds,]
sam = d[!(d$Identifier_1 %in% stds),]

# Calculate the within-run standard variances
vC_rm = vO_rm = n_rm = matrix(NA, length(rf), length(stds))
for(i in seq_along(rf)){
  for(j in seq_along(stds)){
    rm.sub = rm[rm$run == rf[i] & rm$Identifier_1 == stds[j],]
    vC_rm[i, j] = var(rm.sub$d13C.cal)
    vO_rm[i, j] = var(rm.sub$d18O.cal)
    n_rm[i, j] = nrow(rm.sub)
  }
}

# Calculate pooled var for standards
vC_rm = sum((n_rm - 1) * vC_rm, na.rm = TRUE) / (sum(n_rm) - length(stds))
vO_rm = sum((n_rm - 1) * vO_rm, na.rm = TRUE) / (sum(n_rm) - length(stds))

# Calculate the sample variances, each sample replicated in same run
sams = unique(sam$Identifier_1)
vC_sam = vO_sam = n_sam = numeric(length(sams))
for(i in seq_along(sams)){
  sam.sub = sam[sam$Identifier_1 == sams[j],]
  vC_sam[i] = var(sam.sub$d13C.cal)
  vO_sam[i] = var(sam.sub$d18O.cal)
  n_sam[i] = nrow(sam.sub)
}

# Calculate pooled var for samples
vC_sam = sum((n_sam - 1) * vC_sam, na.rm = TRUE) / (sum(n_sam) - length(sams))
vO_sam = sum((n_sam - 1) * vO_sam, na.rm = TRUE) / (sum(n_sam) - length(sams))

# Calculate the random error
muC_rw = sqrt(vC_rm + vC_sam)
muO_rw = sqrt(vO_rm + vO_sam)

# Calculate the MSE for check standard
mar = rm[rm$Identifier_1 == "Marble",]
marC = marO = numeric(length(rf))

for(i in seq_along(rf)){
  mar.sub = mar[mar$run == rf[i],]
  marC[i] = mean(mar.sub$d13C.cal) - 1.9
  marO[i] = mean(mar.sub$d18O.cal) + 11.3
}

mseC = sum(marC ^ 2) / length(marC)
mseO = sum(marO ^ 2) / length(marO)

# Propigate with variance of known value
muC_bias = sqrt(mseC + 0.08 ^ 2)
muO_bias = sqrt(mseO + 0.05 ^ 2)

# Combined uncertainty
muC = sqrt(muC_rw ^ 2 + muC_bias ^ 2)
muO = sqrt(muO_rw ^ 2 + muO_bias ^ 2)
