# Formatting datafile to DMC-appropriate format

#set wd
#setwd("~/Desktop/project/EXG3-prep/controls/wave3/w3_controls_raw_output")
#setwd("~/Desktop/project/EXG3-prep/controls/wave4/w4_controls_raw_output")
#setwd("~/Desktop/project/EXG3-prep/controls/wave5/w5_controls_raw_output")

#load up some libraries:

# Load packages
source('Functions/sourcePkgs.R')

dir.create("~/Desktop/project/EXG3-prep/controls/EXGdata")

# load raw output file
load('w#_controls_RawData_uncleaned.Rdata')

########################################################

# Convert RT and SSD to seconds 
dat <- input %>%
  select(subject, everything()) %>% mutate(
    rt = rt/1000) %>% mutate(
      ssd = ssd/1000)

# Data columns i want
exgdat <- 
  dat[, c("subject",
          "ssd",
          "stimulus",
          "signal",
          "respons",
          "rt")]

# Re-name in line with EXG3 model
names(exgdat) <- c("s", 
                   "SSD", 
                   "S", 
                   "SS", 
                   "R", 
                   "RT")
## --- S ---
# Make the stimulus variable into a factor with 1=square; 2=circle
exgdat$S <- factor(exgdat$S, levels = 1:2, labels = c('square', 'circle'))
# Name the levels of stim factor: s1=square; s2=circle 
levels(exgdat$S) <- c("s1", "s2")

## --- R ---
# Make the response variable into a factor with 0=NR; 1=square; 2=circle
exgdat$R <- factor(exgdat$R, levels = 0:2, labels = c('NR', 'square', 'circle'))
# Name the levels of resp factor: NR=no response; r1=square; r2=circle
levels(exgdat$R) <- c("NR", "r1", "r2")

## --- SS ---
# Make the Signal variable into a factor
exgdat$SS <- factor(exgdat$SS)
# Name the levels of that factor 
levels(exgdat$SS) <- c("GO", "SS")

## --- SSD ---
# Make SSD on go trials infinite
exgdat$SSD[exgdat$SS == "GO"] <- Inf

## -- RT ---
# Go Omissions labelled with NA
exgdat$RT[exgdat$R == "NR"] <- NA

# Order Factor levels
levels(exgdat$S) <- c("s1", "s2")
levels(exgdat$R) <- c("NR", "r1", "r2")
levels(exgdat$SS) <- c("go", "stop")

# Relocate SSD column to the end
exgdat <- exgdat %>% relocate(SSD, .after = RT)

ctl_w#_dat = exgdat

########################################################

save(ctl_w#_dat, file = "~/Desktop/project/EXG3-prep/controls/EXGdata/w#_control_exgdat.Rdata")

#save working environment
dir.create("~/Desktop/project/EXG3-prep/controls/EXGdata/envirs")
save.image( "~/Desktop/project/EXG3-prep/controls/EXGdata/envirs/w#_control_exgdat.Rdata")

########################################################

# Split by ID and save as list

# datlist = split(exgdat, exgdat$s)

# save.image("w#_control_exglist.Rdata")

########################################################
