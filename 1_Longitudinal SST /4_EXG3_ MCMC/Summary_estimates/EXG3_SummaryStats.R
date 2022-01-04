
# EXG3 Parameter summaries 
# This script calculates the median values for each parameter of interest per subject

setwd("~/Desktop/Merv-files/MERV-PROJECT/EXG3_MODELLING/exg3_mcmcSampling")


# First, load all relevant packages in your R console
# Load packages
source('Functions/sourcePkgs.R')

# Load the EXG3 Model
source ("dmc/dmc.R")
load_model ("EXG-SS", "exgSSprobit.R")


# Load all of the datafiles which have been sampled via the "h.RUN.converge" function
exg3_phase3 <- list.files("./exg3_mcmc","*Phase3.RData", full.names="TRUE")
lapply(exg3_phase3,load,.GlobalEnv)


setwd("./exg3_Parameter_estimates")

# Write out the function
pars.estab <- function(hsamples,parameter){
  out <- do.call(rbind, lapply(hsamples, function(x) {as.vector(x$theta[, parameter, ])
  }))
  return(out)
}


## ADHD

# Wave 3
pars <- list()
for (i in colnames(exg3.samples.adhd.w3.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.adhd.w3.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_adhd_wave3_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.adhd.w3.1)
write.csv(SubjectData,file = "adhd_wave3_exg3_params.csv",sep="")



# Wave 4
pars <- list()
for (i in colnames(exg3.samples.adhd.w4.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.adhd.w4.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_adhd_wave4_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.adhd.w4.1)
write.csv(SubjectData,file = "adhd_wave4_exg3_params.csv",sep="")




# Wave 5
pars <- list()
for (i in colnames(exg3.samples.adhd.w5.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.adhd.w5.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_adhd_wave5_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.adhd.w5.1)
write.csv(SubjectData,file = "adhd_wave5_exg3_params.csv",sep="")


#############################################################


# CONTROLS 

# Wave 3
pars <- list()
for (i in colnames(exg3.samples.ctl.w3.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.ctl.w3.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_ctl_wave3_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.ctl.w3.1)
write.csv(SubjectData,file = "ctl_wave3_exg3_params.csv",sep="")



# Wave 4
pars <- list()
for (i in colnames(exg3.samples.ctl.w4.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.ctl.w4.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_ctl_wave4_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.ctl.w4.1)
write.csv(SubjectData,file = "ctl_wave4_exg3_params.csv",sep="")




# Wave 5
pars <- list()
for (i in colnames(exg3.samples.ctl.w5.1[[1]]$theta)){ # same parameter names for each participant, so just used first one
  pars[[i]] <- pars.estab(exg3.samples.ctl.w5.1,i)
}

## Create summary parameters ##
pars[["SSRT"]] <- (pars[["muS"]])+(pars[["tauS"]])
pars[["SSRTvar"]] <- ((pars[["sigmaS"]])^2)+((pars[["tauS"]])^2)
pars[["GoRT"]] <- (pars[["mu.true"]])+(pars[["tau.true"]])
pars[["GoRTvar"]] <- ((pars[["sigma.true"]])^2)+((pars[["tau.true"]])^2)
pars[["GoRTerr"]] <- (pars[["mu.false"]])+(pars[["tau.false"]])
pars[["GoRTvarerr"]] <- ((pars[["sigma.false"]])^2)+((pars[["tau.false"]])^2)
pars[["prob_TF"]] <- pnorm(pars[["tf"]])*100
pars[["prob_GF"]] <- pnorm(pars[["gf"]])*100
save(pars,file = "exg3_ctl_wave5_summaryparams.RData")

# Save to csv
subj.summ <- list()
subj.data <- list()
## save off summary stats ##
for (i in names(pars)){
  subj.summ[[i]] <- lapply(1:length(pars$mu.true[,1]), 
                           function (x) median(pars[[i]][x,]))
  
  subj.data[[i]] <- t(data.frame(subj.summ[[i]]))
}
SubjectData <- data.frame(subj.data)
rownames(SubjectData) <- names(exg3.samples.ctl.w5.1)
write.csv(SubjectData,file = "ctl_wave5_exg3_params.csv",sep="")





