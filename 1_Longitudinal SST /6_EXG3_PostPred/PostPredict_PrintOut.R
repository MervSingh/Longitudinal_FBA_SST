# Fits Posterior Predictives
cwd <- getwd()

# Load packages
source('Functions/sourcePkgs.R')

# Load Posterior Predictives
load("EXG3_PosteriorPredictives.RData")

# CONTROLS ___________________________________________________________________________________________________

# Wave 3

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w3.1,sim=pp.sim4,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("ctl_Wave3_PostPredictive.pdf",width = 15,height = 10)
ggplot.RP.dmc(pp4,include.NR=TRUE)  # Figure 13a
# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp4) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp4,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.ctl.w3.1,sim=pp.sim4,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w3.1,sim=pp.sim4,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(ctl.w3.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.ctl.w3.1[[i]]$data,sim=pp.sim4[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w3.1[[i]]$data,sim=pp.sim4[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp4[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()


#####


# Wave 4

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w4.1,sim=pp.sim5,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("ctl_Wave4_PostPredictive.pdf",width = 15,height = 10)

ggplot.RP.dmc(pp5,include.NR=TRUE)  # Figure 13a

# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp5) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp5,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.ctl.w4.1,sim=pp.sim5,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w4.1,sim=pp.sim5,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(ctl.w4.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.ctl.w4.1[[i]]$data,sim=pp.sim5[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w4.1[[i]]$data,sim=pp.sim5[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp5[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()


#####


# Wave 5

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w5.1,sim=pp.sim6,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("ctl_Wave5_PostPredictive.pdf",width = 15,height = 10)
ggplot.RP.dmc(pp6,include.NR=TRUE)  # Figure 13a
# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp6) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp6,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.ctl.w5.1,sim=pp.sim6,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w5.1,sim=pp.sim6,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(ctl.w5.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.ctl.w5.1[[i]]$data,sim=pp.sim6[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.ctl.w5.1[[i]]$data,sim=pp.sim6[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp6[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()



# ADHD ___________________________________________________________________________________________________

# Wave 3

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w3.1,sim=pp.sim1,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("adhd_Wave3_PostPredictive.pdf",width = 15,height = 10)
ggplot.RP.dmc(pp1,include.NR=TRUE)  # Figure 13a
# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp1) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp1,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.adhd.w3.1,sim=pp.sim1,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w3.1,sim=pp.sim1,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(adhd.w3.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.adhd.w3.1[[i]]$data,sim=pp.sim1[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w3.1[[i]]$data,sim=pp.sim1[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp1[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()


#####


# Wave 4

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w4.1,sim=pp.sim2,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("adhd_Wave4_PostPredictive.pdf",width = 15,height = 10)

ggplot.RP.dmc(pp2,include.NR=TRUE)  # Figure 13a

# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp2) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp2,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.adhd.w4.1,sim=pp.sim2,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w4.1,sim=pp.sim2,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(adhd.w4.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.adhd.w4.1[[i]]$data,sim=pp.sim2[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w4.1[[i]]$data,sim=pp.sim2[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp2[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()


#####


# Wave 5

#set ggplot theme
theme_set(theme_simple())
tmp <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w5.1,sim=pp.sim3,n.intervals=20,
                        violin=.5,percentile.av=FALSE)
tmp <- round(attr(tmp,"cuts"),2)
xlabels <- paste("(",tmp[-length(tmp)],",",tmp[-1],"]",sep="")

pdf("adhd_Wave5_PostPredictive.pdf",width = 15,height = 10)
ggplot.RP.dmc(pp3,include.NR=TRUE)  # Figure 13a
# Next plot the RT distributions, by default as the 10th, 50th and 90th
# percentiles. Note that the NR response level is dropped because by definition
# it has no RT data.
ggplot.RT.dmc(pp3) # Figure 13b
# Plot average cdf (Figure 13c)
plot.pp.dmc(pp3,layout=c(2,2),ylim=c(0,1),fits.pcol="grey",
            model.legend=FALSE,dname="",aname="",x.min.max=c(.25,2))
plot.new()
layout(1)
H_IFTable <- plot_SS_if.dmc(data=exg3.samples.adhd.w5.1,sim=pp.sim3,n.intervals=20,violin=.25,xlab="SSD")
title <- textGrob("Posterior predictive model check for Inhibition Function") ### add the full title from BEESTS
t <- tableGrob(t(H_IFTable),rows = c("SSD (s)","n delays","p value"))
grid.arrange(title,t,ncol=1)

plot.new()
layout(1)
H_ssrtTable <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w5.1,sim=pp.sim3,n.intervals=20,
                                violin=.5,percentile.av=FALSE,xlabels=xlabels,xlab="SSD (s)")
title <- textGrob("Group level posterior predictive model check for Signal Respond RT") ### add the full title from BEESTS
t <- tableGrob(t(H_ssrtTable),rows = c("SSD (s)","n delays","p value","N. Sims"))
grid.arrange(title,t,ncol=1)

# for (i in seq_along(adhd.w5.model2)){
#   plot.new()
#   layout(1)
#   IF_table <- plot_SS_if.dmc(data=exg3.samples.adhd.w5.1[[i]]$data,sim=pp.sim3[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Inhibition Function, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(IF_table[IF_table$n>5,]),rows = c("SSD (s)","n delays","p value"))
#   grid.arrange(title,t,ncol=1)
#   srrt_table <- plot_SS_srrt.dmc(data=exg3.samples.adhd.w5.1[[i]]$data,sim=pp.sim3[[i]],main=paste("Subject",i,sep=" "),minN = 6)
#   title <- textGrob(paste("Posterior predictive model check for Signal Respond RT, Subject",i,sep=" ")) ### add the full title from BEESTS
#   t <- tableGrob(t(srrt_table[srrt_table$nrt>5,])[-4,],rows = c("SSD (s)","n trials","p value"))
#   grid.arrange(title,t,ncol=1)
# 
#   plot.pp.dmc(pp3[[i]],layout=c(2,2),fits.pcol = "red",fits.pch = 20,
#               model.legend=F,dname="",aname="",x.min.max=c(.25,2),pos='topleft',style = 'cdf')
# 
# }
dev.off()