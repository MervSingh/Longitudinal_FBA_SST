#### exg3 ####

# MCMC sampling data examination
# This script relies on the Dynamic Models of Choice (DMC) package in R to interrogate model diagnostics
# and extract summary statistics from the parameters of interest

# First, load all relevant packages in your R console
# Load packages
source('Functions/sourcePkgs.R')

# Load the EXG3 Model
source ("dmc/dmc.R")
load_model ("EXG-SS", "exgSSprobit.R")


# Load all of the datafiles which have been sampled via the "h.RUN.converge" function
exg3_phase3 <- list.files("./exg3_mcmc","*Phase3.RData", full.names="TRUE")
lapply(exg3_phase3,load,.GlobalEnv)


## MODEL SUMMARIES

dir.create("exg3_Parameter_estimates")

# ADHD

# Wave 3
pdf("./exg3_Parameter_estimates/ADHD_Wave3_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.adhd.w3.1)) {
  # Summary stats for each subj
  adhd_w3_summary <- data.frame(summary.dmc(exg3.samples.adhd.w3.1[[i]])[1])
  adhd_w3_quantiles <- data.frame(summary.dmc(exg3.samples.adhd.w3.1[[i]])[2])
  colnames(adhd_w3_summary) <- c("w3_Mean", "w3_SD", "w3_NaiveSD", "w3_Time-seriesSE")
  colnames(adhd_w3_quantiles) <- c("w3_2.5%", "w3_25%", "w3_50%", "w3_75%", "w3_97.5%")
  adhd_w3_df <- cbind(adhd_w3_summary, adhd_w3_quantiles)
  adhd_w3_t <-tableGrob(format(adhd_w3_df, digits = 2))
  adhd_w3_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(adhd_w3_title, adhd_w3_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.adhd.w3.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()


# Wave 4
pdf("./exg3_Parameter_estimates/ADHD_Wave4_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.adhd.w4.1)) {
  # Summary stats for each subj
  adhd_w4_summary <- data.frame(summary.dmc(exg3.samples.adhd.w4.1[[i]])[1])
  adhd_w4_quantiles <- data.frame(summary.dmc(exg3.samples.adhd.w4.1[[i]])[2])
  colnames(adhd_w4_summary) <- c("w4_Mean", "w4_SD", "w4_NaiveSD", "w4_Time-seriesSE")
  colnames(adhd_w4_quantiles) <- c("w4_2.5%", "w4_25%", "w4_50%", "w4_75%", "w4_97.5%")
  adhd_w4_df <- cbind(adhd_w4_summary, adhd_w4_quantiles)
  adhd_w4_t <-tableGrob(format(adhd_w4_df, digits = 2))
  adhd_w4_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(adhd_w4_title, adhd_w4_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.adhd.w4.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()

# Wave 5
pdf("./exg3_Parameter_estimates/ADHD_Wave5_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.adhd.w5.1)) {
  # Summary stats for each subj
  adhd_w5_summary <- data.frame(summary.dmc(exg3.samples.adhd.w5.1[[i]])[1])
  adhd_w5_quantiles <- data.frame(summary.dmc(exg3.samples.adhd.w5.1[[i]])[2])
  colnames(adhd_w5_summary) <- c("w5_Mean", "w5_SD", "w5_NaiveSD", "w5_Time-seriesSE")
  colnames(adhd_w5_quantiles) <- c("w5_2.5%", "w5_25%", "w5_50%", "w5_75%", "w5_97.5%")
  adhd_w5_df <- cbind(adhd_w5_summary, adhd_w5_quantiles)
  adhd_w5_t <-tableGrob(format(adhd_w5_df, digits = 2))
  adhd_w5_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(adhd_w5_title, adhd_w5_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.adhd.w5.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()



#CONTROL


# Wave 3
pdf("./exg3_Parameter_estimates/Control_Wave3_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.ctl.w3.1)) {
  # Summary stats for each subj
  ctl_w3_summary <- data.frame(summary.dmc(exg3.samples.ctl.w3.1[[i]])[1])
  ctl_w3_quantiles <- data.frame(summary.dmc(exg3.samples.ctl.w3.1[[i]])[2])
  colnames(ctl_w3_summary) <- c("w3_Mean", "w3_SD", "w3_NaiveSD", "w3_Time-seriesSE")
  colnames(ctl_w3_quantiles) <- c("w3_2.5%", "w3_25%", "w3_50%", "w3_75%", "w3_97.5%")
  ctl_w3_df <- cbind(ctl_w3_summary, ctl_w3_quantiles)
  ctl_w3_t <-tableGrob(format(ctl_w3_df, digits = 2))
  ctl_w3_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(ctl_w3_title, ctl_w3_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.ctl.w3.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()


# Wave 4
pdf("./exg3_Parameter_estimates/Control_Wave4_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.ctl.w4.1)) {
  # Summary stats for each subj
  ctl_w4_summary <- data.frame(summary.dmc(exg3.samples.ctl.w4.1[[i]])[1])
  ctl_w4_quantiles <- data.frame(summary.dmc(exg3.samples.ctl.w4.1[[i]])[2])
  colnames(ctl_w4_summary) <- c("w4_Mean", "w4_SD", "w4_NaiveSD", "w4_Time-seriesSE")
  colnames(ctl_w4_quantiles) <- c("w4_2.5%", "w4_25%", "w4_50%", "w4_75%", "w4_97.5%")
  ctl_w4_df <- cbind(ctl_w4_summary, ctl_w4_quantiles)
  ctl_w4_t <-tableGrob(format(ctl_w4_df, digits = 2))
  ctl_w4_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(ctl_w4_title, ctl_w4_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.ctl.w4.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()

# Wave 5
pdf("./exg3_Parameter_estimates/Control_Wave5_exg3_Summaries.pdf", width = 11, height = 8)
for (i in seq_along(exg3.samples.ctl.w5.1)) {
  # Summary stats for each subj
  ctl_w5_summary <- data.frame(summary.dmc(exg3.samples.ctl.w5.1[[i]])[1])
  ctl_w5_quantiles <- data.frame(summary.dmc(exg3.samples.ctl.w5.1[[i]])[2])
  colnames(ctl_w5_summary) <- c("w5_Mean", "w5_SD", "w5_NaiveSD", "w5_Time-seriesSE")
  colnames(ctl_w5_quantiles) <- c("w5_2.5%", "w5_25%", "w5_50%", "w5_75%", "w5_97.5%")
  ctl_w5_df <- cbind(ctl_w5_summary, ctl_w5_quantiles)
  ctl_w5_t <-tableGrob(format(ctl_w5_df, digits = 2))
  ctl_w5_title <-
    text_grob(paste("Summary statistics, Subject", i, sep = " "))
  grid.arrange(ctl_w5_title, ctl_w5_t)
  
  # Chains and posteriors for each subject
  plot.dmc(
    exg3.samples.ctl.w5.1[[i]],
    layout = c(3, 2),
    density = T
  )
  plot.new()
  layout(1)
}

dev.off()


save.image("./exg3_Parameter_estimates/exg3_ModelSummaries.RData")
