# Plots SST Performance data - Posterior distributions
cwd <- getwd()
dir.create("./Output")

# Load packages
source('Functions/sourcePkgs.R')


# ADHD ___________________________________________________________________________________________________

# Wave 3

load("/Volumes/Backup Plus/MERV-PROJECT/EXG_7_Plots/samples/exg3_adhd_wave3_summaryparams.RData")

# SSRT #
SSRT <- data.frame(pars$SSRT)
SSRT.long <- melt(t(SSRT))
SSRT.long$Var2 <- as.factor(SSRT.long$Var2)

SSRT.plot <- ggplot(SSRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Stop Signal Reaction Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
SSRT.plot
ggsave('Output/adhd_Wave3_SSRT_subj.jpg')


# GoRT #
GoRT <- data.frame(pars$GoRT)
GoRT.long <- melt(t(GoRT))
GoRT.long$Var2 <- as.factor(GoRT.long$Var2)

GoRT.plot <- ggplot(GoRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Matching Finishing Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
GoRT.plot
ggsave('Output/adhd_Wave3_GoRT_subj.jpg')

# prob_TF #
prob_TF <- data.frame(pars$prob_TF)
prob_TF.long <- melt(t(prob_TF))
prob_TF.long$Var2 <- as.factor(prob_TF.long$Var2)

prob_TF.plot <- ggplot(prob_TF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Trigger Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_TF.plot
ggsave('Output/adhd_Wave3_prob_TF_subj.jpg')

# prob_GF #
prob_GF <- data.frame(pars$prob_GF)
prob_GF.long <- melt(t(prob_GF))
prob_GF.long$Var2 <- as.factor(prob_GF.long$Var2)

prob_GF.plot <- ggplot(prob_GF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Go Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_GF.plot
ggsave('Output/adhd_Wave3_prob_GF_subj.jpg')


# MuS #
MuS <- data.frame(pars$muS)
MuS.long <- melt(t(MuS))
MuS.long$Var2 <- as.factor(MuS.long$Var2)

MuS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Mu for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
MuS.plot
ggsave('Output/adhd_Wave3_MuS_subj.jpg')

# sigmaS #
sigmaS <- data.frame(pars$sigmaS)
sigmaS.long <- melt(t(sigmaS))
sigmaS.long$Var2 <- as.factor(sigmaS.long$Var2)

sigmaS.plot <- ggplot(sigmaS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Sigma for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
sigmaS.plot
ggsave('Output/adhd_Wave3_sigmaS_subj.jpg')


# TauS #
TauS <- data.frame(pars$tauS)
TauS.long <- melt(t(TauS))
TauS.long$Var2 <- as.factor(TauS.long$Var2)

TauS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Tau for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
TauS.plot
ggsave('Output/adhd_Wave3_TauS_subj.jpg')


###################################################################################

rm(list=ls())

# Wave 4

# Load packages
source('Functions/sourcePkgs.R')

load("/Volumes/Backup Plus/MERV-PROJECT/EXG_7_Plots/samples/exg3_adhd_wave4_summaryparams.RData")

# SSRT #
SSRT <- data.frame(pars$SSRT)
SSRT.long <- melt(t(SSRT))
SSRT.long$Var2 <- as.factor(SSRT.long$Var2)

SSRT.plot <- ggplot(SSRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Stop Signal Reaction Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
SSRT.plot
ggsave('Output/adhd_Wave4_SSRT_subj.jpg')


# GoRT #
GoRT <- data.frame(pars$GoRT)
GoRT.long <- melt(t(GoRT))
GoRT.long$Var2 <- as.factor(GoRT.long$Var2)

GoRT.plot <- ggplot(GoRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Matching Finishing Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
GoRT.plot
ggsave('Output/adhd_Wave4_GoRT_subj.jpg')

# prob_TF #
prob_TF <- data.frame(pars$prob_TF)
prob_TF.long <- melt(t(prob_TF))
prob_TF.long$Var2 <- as.factor(prob_TF.long$Var2)

prob_TF.plot <- ggplot(prob_TF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Trigger Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_TF.plot
ggsave('Output/adhd_Wave4_prob_TF_subj.jpg')

# prob_GF #
prob_GF <- data.frame(pars$prob_GF)
prob_GF.long <- melt(t(prob_GF))
prob_GF.long$Var2 <- as.factor(prob_GF.long$Var2)

prob_GF.plot <- ggplot(prob_GF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Go Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_GF.plot
ggsave('Output/adhd_Wave4_prob_GF_subj.jpg')


# MuS #
MuS <- data.frame(pars$muS)
MuS.long <- melt(t(MuS))
MuS.long$Var2 <- as.factor(MuS.long$Var2)

MuS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Mu for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
MuS.plot
ggsave('Output/adhd_Wave4_MuS_subj.jpg')

# sigmaS #
sigmaS <- data.frame(pars$sigmaS)
sigmaS.long <- melt(t(sigmaS))
sigmaS.long$Var2 <- as.factor(sigmaS.long$Var2)

sigmaS.plot <- ggplot(sigmaS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Sigma for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
sigmaS.plot
ggsave('Output/adhd_Wave4_sigmaS_subj.jpg')


# TauS #
TauS <- data.frame(pars$tauS)
TauS.long <- melt(t(TauS))
TauS.long$Var2 <- as.factor(TauS.long$Var2)

TauS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Tau for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
TauS.plot
ggsave('Output/adhd_Wave4_TauS_subj.jpg')


###################################################################################

rm(list=ls())

# Wave 5

# Load packages
source('Functions/sourcePkgs.R')

load("/Volumes/Backup Plus/MERV-PROJECT/EXG_7_Plots/samples/exg3_adhd_wave5_summaryparams.RData")


# SSRT #
SSRT <- data.frame(pars$SSRT)
SSRT.long <- melt(t(SSRT))
SSRT.long$Var2 <- as.factor(SSRT.long$Var2)

SSRT.plot <- ggplot(SSRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Stop Signal Reaction Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
SSRT.plot
ggsave('Output/adhd_Wave5_SSRT_subj.jpg')


# GoRT #
GoRT <- data.frame(pars$GoRT)
GoRT.long <- melt(t(GoRT))
GoRT.long$Var2 <- as.factor(GoRT.long$Var2)

GoRT.plot <- ggplot(GoRT.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Matching Finishing Time",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
GoRT.plot
ggsave('Output/adhd_Wave5_GoRT_subj.jpg')

# prob_TF #
prob_TF <- data.frame(pars$prob_TF)
prob_TF.long <- melt(t(prob_TF))
prob_TF.long$Var2 <- as.factor(prob_TF.long$Var2)

prob_TF.plot <- ggplot(prob_TF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Trigger Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_TF.plot
ggsave('Output/adhd_Wave5_prob_TF_subj.jpg')

# prob_GF #
prob_GF <- data.frame(pars$prob_GF)
prob_GF.long <- melt(t(prob_GF))
prob_GF.long$Var2 <- as.factor(prob_GF.long$Var2)

prob_GF.plot <- ggplot(prob_GF.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Probability Go Failure",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
prob_GF.plot
ggsave('Output/adhd_Wave5_prob_GF_subj.jpg')


# MuS #
MuS <- data.frame(pars$muS)
MuS.long <- melt(t(MuS))
MuS.long$Var2 <- as.factor(MuS.long$Var2)

MuS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Mu for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
MuS.plot
ggsave('Output/adhd_Wave5_MuS_subj.jpg')

# sigmaS #
sigmaS <- data.frame(pars$sigmaS)
sigmaS.long <- melt(t(sigmaS))
sigmaS.long$Var2 <- as.factor(sigmaS.long$Var2)

sigmaS.plot <- ggplot(sigmaS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Sigma for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
sigmaS.plot
ggsave('Output/adhd_Wave5_sigmaS_subj.jpg')


# TauS #
TauS <- data.frame(pars$tauS)
TauS.long <- melt(t(TauS))
TauS.long$Var2 <- as.factor(TauS.long$Var2)

TauS.plot <- ggplot(MuS.long, aes(x=value, colour=Var2)) +
  geom_line(stat='density')+
  theme_classic(base_size = 30,base_family = 'serif')+
  labs(x = "Tau for Stopping Performance",
       y = "Density",colour = "")+
  # scale_colour_hue(h=c(180,181)+15)+
  scale_colour_grey()+
  theme(legend.position="none")
TauS.plot
ggsave('Output/adhd_Wave5_TauS_subj.jpg')

# save workspace
save.image('Output/ADHDPlots.RData')
