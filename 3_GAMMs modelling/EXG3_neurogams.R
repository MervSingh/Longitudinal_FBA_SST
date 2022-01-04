
# GAMS analysis for FBA DATA - TOTAL


#############################################################

# ENVIRONMENT SET-UP

#############################################################

# Load packages
packages <- c("readxl","ggpubr","Hmisc","tidyverse","outliers","ggplot2","parallel","data.table","nlme","mgcv","plyr","dplyr","broom","SemiPar","itsadug","skimr", "fitdistrplus", "mice")
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
dir.create("./FBA_gams_Total")
dir.create("./FBA_gams_Controls_only")


#############################################################

# DATA SET-UP

#############################################################


# Load the QCed FBA subsample for the total sample
total_FBA  <- read_xlsx("total_fba_cleaned.xlsx")

# Format data - ensure all categorical vars (inc "SID) are factors; all continuous vars are numeric
total_FBA <- mutate_at(total_FBA, c(1:6,12,32,52), as.factor) %>%
  mutate_at(c(9:11,13:31,33:51,53), as.numeric)

# make clean dataframe for analysis
dat <- total_FBA

# set NAs
dat <- dat %>% mutate_all(na_if,"NA")

# Check for missingness
sapply(dat, function(x) sum(is.na(x)))

# Create a subset of Data removing all obs with missing etiv
subdat <- subset(dat, (!is.na(dat[["etiv.long"]])))
sapply(subdat, function(x) sum(is.na(x)))

# SIDs that were missing etiv
anti_join(dat, subdat, by= c("SID", "Wave"))

# Filter out ADHD subjects: FOR STUDY 2
condat <- subset(subdat, subdat$group=="CONTROL")
# drop unused "ADHD" factor levels from the CONTROL data
condat$group <- levels(droplevels(condat$group))
condat$group <- as.factor(condat$group)

# Summary stats by Wave x group: TOTAL DATA
neurototal_summ <- subdat %>% dplyr::group_by(Wave) %>% dplyr::group_by(group, .add = T) %>% skimr::skim_without_charts()

# Summary stats by Wave x group: CONTROL DATA
neurocon_summ <- condat %>% dplyr::group_by(Wave) %>% dplyr::group_by(group, .add = T) %>% skimr::skim_without_charts()

# how many unique IDs completed each wave?
ids_con <- subdat %>% filter(group=="CONTROL") %>% group_by(SID) %>% dplyr::summarise(N = n()) %>% mutate(N = as.factor(N))
skim(ids_con)

ids_adhd <- subdat %>% filter(group=="ADHD") %>% group_by(SID) %>% dplyr::summarise(N = n()) %>% mutate(N = as.factor(N))
skim(ids_adhd)

# Write to csv

# TOTAL DATA
write_csv(subdat, "./FBA_gams_Total/FBA_total_final.csv")
write_csv(neurototal_summ, "./FBA_gams_Total/FBA_subsample_TotalDescriptives.csv")
# CONTROL DATA
write_csv(condat, "./FBA_gams_Controls_only/FBA_control_final.csv")
write_csv(neurocon_summ, "./FBA_gams_Controls_only/FBA_subsample_ControlsDescriptives.csv")


# Split dataframe into LHem vs RHem & remove SST vars + b2800 vars

# TOTAL DATA
RHem <- subdat[-c(13:33,34:41)]
LHem <- subdat[-c(13:32,42:50)]

# CONTROL DATA
RHem_con <- condat[-c(13:33,34:41)]
LHem_con <- condat[-c(13:32,42:50)]


# mean centre age, SES, eTIV and meanFWD

# TOTAL DATA
RHem_c <- RHem %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>%
  mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES) %>%
  mutate(etiv.long_c = etiv.long - mean(etiv.long)) %>% relocate(etiv.long_c , .after = etiv.long) %>%
  mutate(meanFWD_c = meanFWD - mean(meanFWD)) %>% relocate(meanFWD_c , .after = meanFWD)

LHem_c <- LHem %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>%
  mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES) %>%
  mutate(etiv.long_c = etiv.long - mean(etiv.long)) %>% relocate(etiv.long_c , .after = etiv.long) %>%
  mutate(meanFWD_c = meanFWD - mean(meanFWD)) %>% relocate(meanFWD_c , .after = meanFWD)

# CONTROL DATA
RHemCON_c <- RHem_con %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>%
  mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES) %>%
  mutate(etiv.long_c = etiv.long - mean(etiv.long)) %>% relocate(etiv.long_c , .after = etiv.long) %>%
  mutate(meanFWD_c = meanFWD - mean(meanFWD)) %>% relocate(meanFWD_c , .after = meanFWD)

LHemCON_c <- LHem_con %>% mutate(age_c = age - mean(age)) %>% relocate(age_c, .after = age) %>%
  mutate(SES_c = SES - mean(SES)) %>% relocate(SES_c, .after = SES) %>%
  mutate(etiv.long_c = etiv.long - mean(etiv.long)) %>% relocate(etiv.long_c , .after = etiv.long) %>%
  mutate(meanFWD_c = meanFWD - mean(meanFWD)) %>% relocate(meanFWD_c , .after = meanFWD)

#############################################################

# GAMS MODELLING

#############################################################

### NEUROIMAGING ###

## Set up for looping

# Create list of vars for lapply to loop different GAM models through

# TOTAL DATA
lh_dvs <- colnames(LHem_c[c(13:18)])
rh_dvs  <- colnames(RHem_c[c(13:18)])

# CONTROL DATA
lh_dvs_CON <- colnames(LHemCON_c[c(13:18)])
rh_dvs_CON <- colnames(RHemCON_c[c(13:18)])


## Transform data to extra long format for lapply to loop over

# TOTAL DATA
lh_xlong <- LHem_c %>% gather(lh_dvs, value, -SID, -SES, -Wave, -wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES_c, -scanner, -meanFWD, -meanFWD_c, -etiv.long, -etiv.long_c)
rh_xlong <- RHem_c %>% gather(rh_dvs, value, -SID, -SES, -Wave, -wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES_c, -scanner, -meanFWD, -meanFWD_c, -etiv.long, -etiv.long_c)

# write to csv
write_csv(lh_xlong, "./FBA_gams_Total/FBA_TotalXlong_LH.csv")
write_csv(rh_xlong, "./FBA_gams_Total/FBA_TotalXlong_RH.csv")


# CONTROL DATA
lh_xlong_CON <- LHemCON_c %>% gather(lh_dvs_CON, value, -SID, -SES, -Wave, -wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES_c, -scanner, -meanFWD, -meanFWD_c, -etiv.long, -etiv.long_c)
rh_xlong_CON <- RHemCON_c %>% gather(rh_dvs_CON, value, -SID, -SES, -Wave, -wave, -sex, -hand, -group, -med, -age, -age_c, -IQ, -SES_c, -scanner, -meanFWD, -meanFWD_c, -etiv.long, -etiv.long_c)

# write to csv
write_csv(lh_xlong_CON, "./FBA_gams_Controls_only/FBA_ConXlong_LH.csv")
write_csv(rh_xlong_CON, "./FBA_gams_Controls_only/FBA_ConXlong_RH.csv")


# Run the GAMS modelling code for all DVs via sourcing the script

# TOTAL DATA
source("./FBA_Source_scripts/LH_fba_GAMSscript_total.R")
source("./FBA_Source_scripts/RH_fba_GAMSscript_total.R")

# CONTROL DATA
source("./FBA_Source_scripts/LH_fba_GAMSscript_controls.R")
source("./FBA_Source_scripts/RH_fba_GAMSscript_controls.R")


# SAVE WORKSPACE
save.image("FBA_GAMS.RData")
