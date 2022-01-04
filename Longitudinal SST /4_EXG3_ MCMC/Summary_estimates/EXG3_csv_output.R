# EXG3

# The script creates a master datafile for ADHD and CONTROL group, containing parameter estimates from all 3 waves

setwd("~/Desktop/EXG3_mcmcSampling")

###LOAD PACKAGES AND SET DIRECTORIES
# Load packages
source('Functions/sourcePkgs.R')

setwd("~/Desktop/exg3_mcmcSampling/exg3_Parameter_estimates")

# Load csv files of the summary estimates for ADHD and controls at each wave & rename first column to subject
#ADHD
adhd_w3 = read_csv("adhd_wave3_exg3_params.csv")
adhd_w3 = adhd_w3 %>% rename_with( ~ paste("w3", .x, sep = "_"))
adhd_w3 = adhd_w3 %>% rename(subject = w3_X1)
  

adhd_w4 = read_csv("adhd_wave4_exg3_params.csv")
adhd_w4 = adhd_w4 %>% rename_with( ~ paste("w4", .x, sep = "_"))
adhd_w4 = adhd_w4 %>% rename(subject = w4_X1)


adhd_w5 = read_csv("adhd_wave5_exg3_params.csv")
adhd_w5 = adhd_w5 %>% rename_with( ~ paste("w5", .x, sep = "_"))
adhd_w5 = adhd_w5 %>% rename(subject = w5_X1)

# Combine all ADHD dataframes into a single masterfile
adhd_tmp <- merge(adhd_w3, adhd_w4, all=TRUE)
total_adhd <- merge(adhd_tmp, adhd_w5, all=TRUE)

# Write out csv file
write.csv(total_adhd, file = "ADHD_Total_EXG3.csv")

#########################################################

#CONTROL
ctl_w3 = read_csv("ctl_wave3_exg3_params.csv")
ctl_w3 = ctl_w3 %>% rename_with( ~ paste("w3", .x, sep = "_"))
ctl_w3 = ctl_w3 %>% rename(subject = w3_X1)


ctl_w4 = read_csv("ctl_wave4_exg3_params.csv")
ctl_w4 = ctl_w4 %>% rename_with( ~ paste("w4", .x, sep = "_"))
ctl_w4 = ctl_w4 %>% rename(subject = w4_X1)


ctl_w5 = read_csv("ctl_wave5_exg3_params.csv")
ctl_w5 = ctl_w5 %>% rename_with( ~ paste("w5", .x, sep = "_"))
ctl_w5 = ctl_w5 %>% rename(subject = w5_X1)


# Combine all CONTROL dataframes into a single masterfile
ctl_tmp <- merge(ctl_w3, ctl_w4, all=TRUE)
total_ctl <- merge(ctl_tmp, ctl_w5, all=TRUE)

# Write out csv file
write.csv(total_ctl, file = "CONTROL_Total_EXG3.csv")
