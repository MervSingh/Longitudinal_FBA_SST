# Short script to take posterior predictive samples.

cwd <- getwd()

# Load packages
source('Functions/sourcePkgs.R')

# Load all of the datafiles which have been sampled via the "h.RUN.converge" function
exg3_phase3 <- list.files(path = "./Phase3_samples", pattern = "*Phase3.RData", full.names="TRUE")
lapply(exg3_phase3,load,.GlobalEnv)

dir.create("./PP_output")

### STANDARD GOODNESS OF FIT CHECKS ----

# # Create posterior predictives to check goodness of fit;
# # Generates 1000 simulated datasets based on the assumptions that the fitted model is a true representation of SST behaviour
# # This can be SLOW unless you have lots of cores available

setwd("./PP_output")

# ADHD

# Wave 3
pp1 <- h.post.predict.dmc(exg3.samples.adhd.w3.1, cores = cores, gglist = TRUE, n.post = 1000) # errored out
# Save Data
save(pp1, file = paste(cwd, "ADHD_w3_PP.RData", sep = ""))

# Wave 4
pp2 <- h.post.predict.dmc(exg3.samples.adhd.w4.1, cores = cores, gglist = TRUE, n.post = 1000)
# Save Data
save(pp2, file = paste(cwd, "ADHD_w4_PP.RData", sep = ""))

# Wave 5
pp3 <- h.post.predict.dmc(exg3.samples.adhd.w5.1, cores = cores, gglist = TRUE, n.post = 1000)
# Save Data
save(pp3, file = paste(cwd, "ADHD_w5_PP.RData", sep = ""))


# CONTROL

# Wave 3
pp4 <- h.post.predict.dmc(exg3.samples.ctl.w3.1, cores = cores, gglist = TRUE, n.post = 1000) # errored out
# Save Data
save(pp4, file = paste(cwd, "ctl_w3_PP.RData", sep = ""))

# Wave 4
pp5 <- h.post.predict.dmc(exg3.samples.ctl.w4.1, cores = cores, gglist = TRUE, n.post = 1000)
# Save Data
save(pp5, file = paste(cwd, "ctl_w4_PP.RData", sep = ""))

# Wave 5
pp6 <- h.post.predict.dmc(exg3.samples.ctl.w5.1, cores = cores, gglist = TRUE, n.post = 1000)
# Save Data
save(pp6, file = paste(cwd, "ctl_w5_PP.RData", sep = ""))


setwd(cwd)


#########


# The following plot functions require saving posterior predictve simulations;
# Again this can be SLOW

dir.create("./PP_sim")
setwd("./PP_sim")

# ADHD

# Wave 3
pp.sim1 <- h.post.predict.dmc(exg3.samples.adhd.w3.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim1, file = paste(cwd, "ADHD_w3_PPsim.RData", sep = ""))

# Wave 4
pp.sim2 <- h.post.predict.dmc(exg3.samples.adhd.w4.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim2, file = paste(cwd, "ADHD_w4_PPsim.RData", sep = ""))

# Wave 5
pp.sim3 <- h.post.predict.dmc(exg3.samples.adhd.w5.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim3, file = paste(cwd, "ADHD_w5_PPsim.RData", sep = ""))


# CONTROL

# Wave 3
pp.sim4 <- h.post.predict.dmc(exg3.samples.ctl.w3.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim4, file = paste(cwd, "ctl_w3_PPsim.RData", sep = ""))

# Wave 4
pp.sim5 <- h.post.predict.dmc(exg3.samples.ctl.w4.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim5, file = paste(cwd, "ctl_w4_PPsim.RData", sep = ""))

# Wave 5
pp.sim6 <- h.post.predict.dmc(exg3.samples.ctl.w5.1, cores = cores, save.simulation = TRUE, n.post = 1000)
# Save Data
save(pp.sim6, file = paste(cwd, "ctl_w5_PPsim.RData", sep = ""))


setwd(cwd)

# Save complete file 
save.image("EXG3_PosteriorPredictives.RData")
