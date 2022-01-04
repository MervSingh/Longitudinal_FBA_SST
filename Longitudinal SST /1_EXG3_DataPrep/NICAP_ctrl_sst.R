#set wd
setwd("~/Desktop/project/EXG3-prep/controls/wave3")
#setwd("~/Desktop/project/EXG3-prep/controls/wave4")
#setwd("~/Desktop/project/EXG3-prep/controls/wave5")

#load up some libraries:

# Load packages
source('Functions/sourcePkgs.R')

# Splitting a data frame by a variable and saving separate csv files
# Uncomment the following lines if your raw SST data is saved in a single csv file
# Otherwise, skip to load raw SST files

#list.files() #list the files
#dat <- read.csv(".csv") #load the raw SST csv file - named according to wave
#View(dat) #view the dataframe

#dat2 <- dat #make a copy of the dataframe to work from (I typically do this so as the orignial dataframe is mantained)

#spt2 <- split(dat2, dat2$ID) #create a 3rd dataframe by splitting the 2nd dataframe by ID

#lapply(names(spt2), function(x){
  #write_csv(spt2[[x]], path = paste("stop_", x, ".csv", sep = ""))
#}) #Write each split dataframe into its own csv file and save it into your working directory

##############################

rm(list=ls()) #refreshes the working environment

#Load the raw STOP-IT datafiles into your working directory

input <- data.frame() #create empty dataframe
oldWD <- getwd() #places your current working directory in an object

files <- dir(pattern = "stop_*") #list the raw datafiles
for (i in files) {
    tmp <- read.csv (i, header = T) #read the file into R
  
  #extract subject number
  labels <-  unlist(strsplit(files[match(i, files)], "_"))
  subject_label <- unlist(strsplit(labels[2], "\\."))[1]
  tmp$subject <- as.numeric(subject_label)
  
  #add the content to a single data frame
  input <- rbind (input, tmp)
  rm(tmp)
}

setwd(oldWD) #to go back to main folder


##############################


# Data preparation

#This prepares the raw datafiles and creates new variables which will serve in the calculation of SSRT and other performance measures.

# Remove redundant ID column
input <- input[-c(1)]

# convert the signal variable to a categorical factor
# (0 = 'nosignal' [i.e. go trial]; 1 = 'signal' [i.e. stop trial])
input$signal <- factor(input$signal, levels = 0:1, labels = c('nosignal', 'signal'))

# optional: convert the stim variable to a factor
# (1 = square; 2 = circle)
#input$stim <- factor(input$stim, levels = 1:2, labels = c('square', 'circle'))

# check frequency counts for each catergorical variable:
table(input$subject, input$signal) #no. of nosignal & signal trials per subject 

# 144 go trials & 48 stop-signal trials per subject
# if subject does not have appropriate numbers as above, remove
# input <- subset(input, input$subject != ##)

table(input$subject, input$stim) #no. of  stimuli per subject
table(input$signal, input$stim) #no. of stimuli for each trial-type

# create new variables for calculation of p(correct)
input$acc <- ifelse(input$correct == 2, 1, 0) #accuracy
input$miss <- ifelse(input$respons == 0, 1, 0) #missed
input$resp <- ifelse(input$respons > 0, 1, 0)  #respond

# Design checks for the new variables
table(input$acc, input$correct)
table(input$respons, input$resp)
table(input$respons, input$miss)


##############################


# Calculation of Stop-signal performance

# in the original version, the first trial of a block was excluded ---
# uncomment the next line if you would like to do this again
# input <- subset(input, trial != '0')

# function to analyse stop-signal performance simultaneously for all subjects:
funcSignal <- function(data){
  
  # signal data: prespond, ssd, true SSRT, and signal-respond RT
  signal <- subset(data, signal == 'signal')
  presp <-  mean(signal$resp)
  ssd <- mean(signal$ssd)
  
  signal.resp <- subset(signal, resp == '1')
  signal.resp.rt <- mean(signal.resp$rt)
  
  # subset no signal data: with and without missed responses
  # for the missed responses, set RT to max RT of the subject/condition
  nosignal <- subset(data, signal == 'nosignal')
  pmiss <- 1 - mean(nosignal$resp) # determine the actual probability of a missed go response
  nosignal_resp <- subset(nosignal, resp == '1')
  nosignal$rt <- ifelse(nosignal$rt == 0, max(nosignal_resp$rt), nosignal$rt)
  
  # --- estimate 1 --- all no-signal trials are INcluded when the nth RT is determined
  ## determine nth RT
  nthRT <- quantile(nosignal$rt, probs = presp, type = 6)
  
  ## SSRT(integration) = nthRT - ssd
  SSRTint <- nthRT - ssd
  
  #  --- estimate 2 ---  estimate SSRT with the mean method
  # DO NOT USE; included only for comparison purposes
  mRT <- mean(nosignal_resp$rt)
  SSRTmean <- mRT - ssd
  
  # Also calculate no-signal RT for go trials with a response,
  # and the difference with signal-respond RT
  nosignal.resp.rt <- mean(nosignal_resp$rt)
  race.check <- nosignal.resp.rt - signal.resp.rt
  
  # Return all data -- replace #w# with wave number i.e w1
  return(data.frame(#w#_presp = round(presp,5),
                    #w#_ssd = round(ssd,5),
                    #w#_nthRT = round(nthRT,5),
                    #w#_SSRTint = round(SSRTint,5),
                    #w#_SSRTmean = round(SSRTmean,5),
                    #w#_usRT = round(signal.resp.rt,5),
                    #w#_goRT_all = round(nosignal.resp.rt,5),
                    #w#_raceCheck = round(race.check,5),
                    #w#_goPmiss = round(pmiss,5)))
}


signal.cast <- ddply(input, .(subject), funcSignal) #Analyse SST performance for all subjects using the function


# --- some extra no-signal data ---

# subset data
nosignal.input <- subset(input, signal == 'nosignal')

# create molten object
nosignal.molten <- melt(nosignal.input, id.var = c('subject', 'correct', 'resp'), measure.var = c('acc', 'miss', 'rt'))

# calculate percent correct
# accuracy or p(correct) = correct trials / (correct trials + incorrect trials).
# trials without a response (or anticpatory responses) are omitted.
acc.cast <- cast (nosignal.molten, subject ~ ., mean, subset = variable == "acc" &  resp == "1")
names(acc.cast)[2] <- "acc"

# calculate RT for correct responses
rt.cast <- cast (nosignal.molten,  subject ~ ., mean, subset = variable == "rt" &  correct == "2")
names(rt.cast)[2] <- "rt"

sd.rt.cast <- cast (nosignal.molten,  subject ~ ., sd, subset = variable == "rt" &  correct == "2")
names(sd.rt.cast)[2] <- "sdrt"


# Combine all data and write the analysed output to csv and Rdata files 
# replace #w# with wave number i.e w1
combined <- signal.cast
combined$#w#_goERR <- round(1-acc.cast$acc,5)
combined$#w#_goRT_correct <- round(rt.cast$rt,5)
combined$#w#_goRT_correct_sd <- round(sd.rt.cast$sdrt,5)

dir.create("./w#_controls_processed_data") #create a new directory to place all ouput files
save(combined, file = "./w#_controls_processed_data/w#_controls_Results_uncleaned.Rdata") #Rdata
write.csv(combined, "./w#_controls_processed_data/w#_controls_Results_uncleaned.csv", row.names = F) #csv

#Save the raw input datafile:
dir.create("./w#_controls_raw_output")
save(input, file = "./w#_controls_raw_output/w#_controls_RawData_uncleaned.RData")

save.image("wave#_controls.Rdata")


##############################


# Combine SST data from all waves to a single file


# Read in SST data for each wave
wave3 <- read.csv("./wave#/w#_controls_processed_data/w#_controls_Results_uncleaned.csv")
wave4 <- read.csv("./wave#/w#_controls_processed_data/w#_controls_Results_uncleaned.csv")
wave5 <- read.csv("./wave#/w#_controls_processed_data/w#_controls_Results_uncleaned.csv")

# Merge all SST data into one dataframe
dat <- merge(wave#, wave#, by = "subject"  , all = T)
dat2 <- merge(dat, wave#, by = "subject"  , all = T)

# Save merged SST data
write_csv(dat2, "controls_sst_all.csv")


##############################
