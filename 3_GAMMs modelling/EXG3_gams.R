
# GAMS analysis for EXG3 DATA - FOR STUDIES 2 AND 3 COMBINED

#############################################################

# ENVIRONMENT SET-UP

#############################################################

# setwd("~/Desktop/EXG_6_Stats/GAMS_Modelling/EXG3")
# setwd("~/Desktop/Merv-files/EXTRA-FOLDERS/workingfile-FOR_GAMS/EXG3")

# Load packages
packages <- c("tidyverse","outliers","ggplot2","parallel","data.table","nlme","mgcv","plyr","dplyr","broom","SemiPar","itsadug","skimr", "fitdistrplus")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

# set number of cores
# number of cores selelcted will be equal to the total number of cores on your machine
numcores <- detectCores()

# set select function from dplyr
select <- dplyr::select

# set directory paths for the input and output data
dir.create("./Controls_only")
dir.create("./Total_sample")


#############################################################

# DATA SET-UP

#############################################################

# Load EXG3 masterfile data
Data <- read.csv("../data_setup/outfiles/exg3_Masterfile.csv", header = T)

# Format data - ensure all categorical vars (inc "SID) are factors; all continuous vars are numeric
Data <- mutate_at(Data, vars('SID', 'Wave', 'sex', 'hand', 'med', 'group'), as.factor) %>% 
  mutate_at(vars('age','IQ','SES','mu.true','sigma.true','tau.true','tf',
                 'muS','sigmaS','tauS','gf','SSRT','SSRTvar','GoRT','GoRTvar','prob_TF','prob_GF'), as.numeric)

# Create a subset of Data removing all obs with missing age
Data_subset <- subset(Data, (!is.na(Data[["age"]])))

# Create a subset of Data removing all obs with missing SES
Data_subset2 <- subset(Data_subset, (!is.na(Data_subset[["SES"]])))

# SIDs that were missing ages
anti_join(Data, Data_subset, by= c("SID", "Wave"))

# SIDs that were missing SES
anti_join(Data_subset, Data_subset2, by= c("SID", "Wave"))


#############################################################

# DATA EXPLORATION AND CLEANING

#############################################################

# make clean dataframe for analysis
dat <- Data_subset2

# Create a new variable for the TOTAL dataset by grand mean centering age & SES
Dat_total <- dat %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>% mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES)

# Further subset data to remove all ADHD subjects: FOR STUDY 2
Con <- subset(dat, dat$group=="CONTROL")

# Create a new variable for the CONTROL dataset by grand mean centering age & SES
Con <- Con %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>%
  mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES)


# Convert all SST parameters from secs to msecs

# TOTAL DATA
Dat_total <- Dat_total %>% mutate(muS = muS*1000) %>%
  mutate(sigmaS = sigmaS*1000) %>% 
  mutate(tauS = tauS*1000) %>% 
  mutate(mu.true = mu.true*1000) %>% 
  mutate(sigma.true = sigma.true*1000) %>% 
  mutate(tau.true = tau.true*1000) %>% 
  mutate(mu.false = mu.false*1000) %>%
  mutate(sigma.false = sigma.false*1000) %>% 
  mutate(tau.false = tau.false*1000) %>% 
  mutate(GoRT = GoRT*1000) %>% 
  mutate(GoRTvar = GoRTvar*1000) %>% 
  mutate(GoRTerr = GoRTerr*1000) %>% 
  mutate(GoRTvarerr = GoRTvarerr*1000) %>% 
  mutate(SSRT = SSRT*1000) %>% 
  mutate(SSRTvar = SSRT*1000)

# CONTROL DATA
Con <- Con %>% mutate(muS = muS*1000) %>%
  mutate(sigmaS = sigmaS*1000) %>%
  mutate(tauS = tauS*1000) %>%
  mutate(mu.true = mu.true*1000) %>%
  mutate(sigma.true = sigma.true*1000) %>%
  mutate(tau.true = tau.true*1000) %>%
  mutate(mu.false = mu.false*1000) %>%
  mutate(sigma.false = sigma.false*1000) %>%
  mutate(tau.false = tau.false*1000) %>%
  mutate(GoRT = GoRT*1000) %>%
  mutate(GoRTvar = GoRTvar*1000) %>%
  mutate(GoRTerr = GoRTerr*1000) %>%
  mutate(GoRTvarerr = GoRTvarerr*1000) %>%
  mutate(SSRT = SSRT*1000) %>%
  mutate(SSRTvar = SSRT*1000)

# drop unused "ADHD" factor levels from the CONTROL data
Con$group <- levels(droplevels(Con$group))
Con$group <- as.factor(Con$group)

# Summary stats by Wave x group: TOTAL DATASET
sumStats_total <- Dat_total %>% dplyr::group_by(Wave) %>% dplyr::group_by(group, .add = T) %>% skimr::skim_without_charts()

# Summary stats by Wave x group: CONTROL DATASET
sumStats_con <- Con %>% dplyr::group_by(Wave) %>% dplyr::group_by(group, .add = T) %>% skimr::skim_without_charts()

# how many unique IDs completed each wave?
ids_conn <- Dat_total %>% group_by(SID) %>% filter(group=="CONTROL") %>% dplyr::summarise(N = n()) %>% mutate(N = as.factor(N))
skim(ids_conn)
ids_adhd <- Dat_total %>% group_by(SID) %>% filter(group=="ADHD") %>% dplyr::summarise(N = n()) %>% mutate(N = as.factor(N))
skim(ids_adhd)

# Write to csv
# TOTAL DATA
write_csv(sumStats_total, "./Total_sample/exg3_CombinedDescriptives.csv")
# Save longform data for FBA steps
write_csv(Dat_total, "./Total_sample/exg3_total_longform.csv")


# CONTROL DATA
# Write to csv
write_csv(sumStats_con, "./Controls_only/exg3_ControlsDescriptives.csv")
# Write to csv
write.csv(Con,"./Controls_only/exg3_Controls_final.csv", row.names = F)



#############################################################

# GAMS MODELLING

#############################################################

## Set up for looping

# Create list of vars for lapply to loop different GAM models through

# TOTAL DATA
dvstotal <- colnames(Dat_total[c(12:17,19:21,23,25,29:30)])
# CONTROL DATA
dvscon <- colnames(Con[c(19:21,23)])

## Transform data to extra long format for lapply to loop over

# TOTAL DATA
Dat_xlong <- Dat_total %>% gather(dvstotal, value, -SID, -SES, -Wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES, -SES_c)
# write to csv
write_csv(Dat_xlong, "./Total_sample/exg3_TotalXlong.csv")

# CONTROL DATA
Con_xlong <- Con %>% gather(dvscon, value, -SID, -SES, -Wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES_c)
# write to csv
write_csv(Con_xlong, "./Controls_only/exg3_ConXlong.csv")


# Run the GAMS modelling code for all DVs via sourcing the script
source("./EXG3_Source_scripts/exg3_GAMSscript_total.R")
source("./EXG3_Source_scripts/exg3_GAMSscript_controls.R")


# SAVE WORKSPACE
save.image("GAMS_EXG3_masterfile.RData")
