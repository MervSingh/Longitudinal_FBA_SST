# Figures for STUDY 2 Manuscript


#############################################################

# ENVIRONMENT SET-UP

#############################################################

# Load packages
packages <- c("ggpubr","Hmisc","tidyverse","outliers","ggplot2","parallel","data.table","nlme","mgcv","plyr","dplyr","broom","SemiPar","itsadug","skimr", "fitdistrplus")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)



################## CONTROL DATASET FIGURES #######################################

# Read behavioural and FBA subsampe data into R
behavdat = read_csv("./behavioural-dat/exg3_Controls_final.csv")
neurodat = readxl::read_xlsx("./fbasample-dat/total_fba_cleaned.xlsx")


# set NAs
behavdat <- behavdat %>% mutate_all(na_if,"NA")
neurodat <- neurodat %>% mutate_all(na_if,"NA")


# Check for missingness
sapply(behavdat, function(x) sum(is.na(x)))
sapply(neurodat, function(x) sum(is.na(x)))


# Create a subset of of  neurodat removing all obs with missing etiv - as these were emoved in the FBA GAMS
neurodat <- subset(neurodat, (!is.na(neurodat[["etiv.long"]])))

# Format data - ensure all categorical vars (inc "SID) are factors; all continuous vars are numeric
behavdat <- mutate_at(behavdat, c(1:5,11), as.factor) %>% mutate_at(vars(6:10,12:30), as.numeric)
neurodat <- mutate_at(neurodat, c(1:6,12,32,52), as.factor) %>% mutate_at(c(9:11,13:31,34:51,53), as.numeric)

#Filter out ADHD
neurodat <- subset(neurodat, neurodat$group=="CONTROL")
neurodat$group <- levels(droplevels(neurodat$group))
neurodat$group <- as.factor(neurodat$group)


#################################################################################

# BEHAVIOURAL SAMPLE

################################################################################

# Check distributions for each DV

#STOPPING PARAMTERS
bp1 <- ggplot(behavdat, aes(x=Wave, y=muS, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp2 <- ggplot(behavdat, aes(x=Wave, y=sigmaS, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp3 <- ggplot(behavdat, aes(x=Wave, y=tauS, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp4 <- ggplot(behavdat, aes(x=Wave, y=SSRT, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp5 <- ggplot(behavdat, aes(x=Wave, y=SSRTvar, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)

#TRIGGER & GO FAILURES
bp6 <- ggplot(behavdat, aes(x=Wave, y=prob_TF, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp7 <- ggplot(behavdat, aes(x=Wave, y=prob_GF, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)

# GO MATCHING PARAMETERS
bp8 <- ggplot(behavdat, aes(x=Wave, y=GoRT, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp9 <- ggplot(behavdat, aes(x=Wave, y=,GoRTvar, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp10 <- ggplot(behavdat, aes(x=Wave, y=,mu.true, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp11 <- ggplot(behavdat, aes(x=Wave, y=sigma.true, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp12 <- ggplot(behavdat, aes(x=Wave, y=tau.true, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)

# GO MIS-MATCHING PARAMETERS
bp13 <- ggplot(behavdat, aes(x=Wave, y=,mu.false, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp14 <- ggplot(behavdat, aes(x=Wave, y=sigma.false, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp15 <- ggplot(behavdat, aes(x=Wave, y=tau.false, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp16 <- ggplot(behavdat, aes(x=Wave, y=,GoRTerr, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)
bp17 <- ggplot(behavdat, aes(x=Wave, y=,GoRTvarerr, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3)


# Arrange all images
stoplot = ggarrange(bp4,bp5,bp1,bp2,bp3,bp6,bp7,labels = c("A", "B", "C", "D", "E", "F", "G"),common.legend = TRUE, legend = "bottom")
goplot = ggarrange(bp8,bp9,bp16,bp17,bp10,bp11,bp12,bp13,bp14,bp15,labels = c("A", "B", "C", "D","E","F","G","H","I","J"),common.legend = TRUE, legend = "bottom")

# save images
ggsave("Controls_stop_violplots.png", stoplot)
ggsave("Controls_go_violplots.png", goplot)


#################################################################################

# NEUROIMAGING SAMPLE

################################################################################


# Check distributions for each DV

# LEFT IFG-PRESMA
np1 <- ggplot(neurodat, aes(x=Wave, y=FD_LH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np2 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_LH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np3 <- ggplot(neurodat, aes(x=Wave, y=FDC_LH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

IFG_PRESMA_L = ggarrange(np1,np2,np3,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
IFG_PRESMA_L = annotate_figure(IFG_PRESMA_L, top = text_grob("Left IFG-PRESMA pathway",color = "red", face = "bold", size = 14))
ggsave("IFG_PRESMA_L.png", IFG_PRESMA_L)
dev.off()

# RIGHT IFG-PRESMA
np4 <- ggplot(neurodat, aes(x=Wave, y=FD_RH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np5 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_RH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np6 <- ggplot(neurodat, aes(x=Wave, y=FDC_RH_IFG_PRESMA, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

IFG_PRESMA_R = ggarrange(np4,np5,np6,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
IFG_PRESMA_R = annotate_figure(IFG_PRESMA_R, top = text_grob("Right IFG-PRESMA pathway",color = "red", face = "bold", size = 14))
ggsave("IFG_PRESMA_R.png", IFG_PRESMA_R)
dev.off()



# LEFT IFG-STN
np7 <- ggplot(neurodat, aes(x=Wave, y=FD_LH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np8 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_LH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np9 <- ggplot(neurodat, aes(x=Wave, y=FDC_LH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

IFG_STN_L = ggarrange(np7,np8,np9,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
IFG_STN_L = annotate_figure(IFG_STN_L, top = text_grob("Left IFG-STN pathway",color = "red", face = "bold", size = 14))
ggsave("IFG_STN_L.png", IFG_STN_L)
dev.off()

# RIGHT IFG-STN
np10 <- ggplot(neurodat, aes(x=Wave, y=FD_RH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np11 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_RH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np12 <- ggplot(neurodat, aes(x=Wave, y=FDC_RH_IFG_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

IFG_STN_R = ggarrange(np10,np11,np12,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
IFG_STN_R = annotate_figure(IFG_STN_R, top = text_grob("Right IFG-STN pathway",color = "red", face = "bold", size = 14))
ggsave("IFG_STN_R.png", IFG_STN_R)
dev.off()


#LEFT PRESMA-STN
np13 <- ggplot(neurodat, aes(x=Wave, y=FD_LH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np14 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_LH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np15 <- ggplot(neurodat, aes(x=Wave, y=FDC_LH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

PRESMA_STN_L = ggarrange(np13,np14,np15,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
PRESMA_STN_L = annotate_figure(PRESMA_STN_L, top = text_grob("Left PRESMA-STN pathway",color = "red", face = "bold", size = 14))
ggsave("PRESMA_STN_L.png", PRESMA_STN_L)
dev.off()

#RIGHT PRESMA-STN
np16 <- ggplot(neurodat, aes(x=Wave, y=FD_RH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FD")
np17 <- ggplot(neurodat, aes(x=Wave, y=LOGFC_RH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean logFC")
np18 <- ggplot(neurodat, aes(x=Wave, y=FDC_RH_PRESMA_STN, fill=Wave)) + geom_violin(trim=FALSE, alpha=0.3) + geom_jitter(shape=16, position=position_jitter(0.2), alpha=4/10) + theme_classic() + geom_boxplot(width=0.1, fill="white", alpha=0.3) + ylab("mean FDC")

PRESMA_STN_R = ggarrange(np16,np17,np18,labels = c("A", "B", "C"),common.legend = TRUE,legend = "bottom")
PRESMA_STN_R = annotate_figure(PRESMA_STN_R, top = text_grob("Right PRESMA-STN pathway",color = "red", face = "bold", size = 14))
ggsave("PRESMA_STN_R.png", PRESMA_STN_R)
dev.off()


#################################################################################

# BEHAV + NEURO SAMPLE SPREAD OF SCORES

################################################################################

# order age variable 

# behav data

bplot <- behavdat %>% mutate(age = round(as.numeric(age),2)) %>% arrange(age) %>%
  mutate(SID = factor(SID, unique(SID)))

# neuro data

nplot <- neurodat %>% mutate(age = round(as.numeric(age),2)) %>% arrange(age) %>%
  mutate(SID = factor(SID, unique(SID)))


# plot spread of data

# behav data

plot1 <- ggplot(bplot, aes(y=SID, x=age, group=SID, colour=sex))+
geom_line(size=.6,alpha=0.2) +
ylab("Participants") +  #Specify titles for y-axis...
xlab("Age") +           #x-axis...
geom_point(size=2) +
theme_bw() +
theme_minimal() +
theme(axis.line = element_line(colour = "black"),
axis.text.y = element_blank(),
axis.ticks.y = element_blank(),
legend.position="none",
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
panel.background = element_blank()) + theme(legend.position="right")


# neuro data

plot2 <- ggplot(nplot, aes(y=SID, x=age, group=SID, colour=sex))+
  geom_line(size=.6,alpha=0.2) +
  ylab("Participants") +  #Specify titles for y-axis...
  xlab("Age") +           #x-axis...
  geom_point(size=2) +
  theme_bw() +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(legend.position="right")


# combine plots to one figure 
fig1 = ggarrange(plot1,plot2, labels = c("a", "b"),common.legend = TRUE, legend = "bottom")

# save
ggsave("spread_data_controls.png", fig1)
dev.off()



# plots X scanner for neuro data

splot <- neurodat %>%
  arrange(age) %>%
  mutate(SID = factor(SID, unique(SID)))

scancolour <- c("#999999", "#E69F00")

plotscan <- ggplot(splot, aes(y=SID, x=age, group=SID, colour=scanner)) +
  geom_line(size=.6,alpha=0.6) + 
  ylab("Participants") +
  xlab("Age") +
  labs(colour = "Scanner") +
  geom_point(size=2) +
  theme_bw() +
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values=scancolour)

# save
ggsave("spread_data_scan_con.png", plotscan)
dev.off()


######################################################################

# SAVE ENVIRONMENT

######################################################################

save.image("Study2_figures.RData")

######################################################################
