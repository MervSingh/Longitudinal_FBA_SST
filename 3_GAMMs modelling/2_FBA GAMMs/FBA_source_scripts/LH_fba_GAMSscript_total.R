
# Study 3 TOTAL - LEFT HEMISPHERIC FRONTO-BASAL-GANGLIA FBA TRAJECTORIES

dir.create("./FBA_gams_Total/LH")

# MODELS WITH COVARS (FULL MODEL)

# GAM models for EXG3 parameters
dir.create("./FBA_gams_Total/LH/fullmodel")

models<-lapply(X=as.character(as.list(lh_dvs)),
               df=lh_xlong,
               FUN=function(dv_name, df) {
                 print(dv_name)
                 
                 # now create a new dataframe within the loop that is filtered to dv_name
                 adf<-df %>% filter(lh_dvs==dv_name) %>%
                 mutate(SID = as.factor(SID),
                        Wave = as.factor(Wave), 
                        scanner = as.factor(scanner),
                        etiv.long_c = as.numeric(etiv.long_c),
                        meanFWD_c = as.numeric(meanFWD_c),
                        age_c = as.numeric(age_c),
                        sex = as.factor(sex),
                        hand = as.factor(hand), 
                        med = as.factor(med), 
                        SES_c = as.numeric(SES_c),
                        group = as.factor(group))
  
                 # ordered factoring of sex
                 adf$OFsex <- as.factor(adf$sex)
                 # change factor to ordered factor:
                 adf$OFsex <- as.ordered(adf$OFsex)
                 # change contrast to treatment coding (difference curves)
                 contrasts(adf$OFsex) <- 'contr.treatment'
                 # Inspect contrasts:
                 contrasts(adf$OFsex)
                 
        
                 # ordered factoring of scanner type
                 adf$OFscanner <- as.factor(adf$scanner)
                 # change factor to ordered factor:
                 adf$OFscanner <- as.ordered(adf$OFscanner)
                 # change contrast to treatment coding (difference curves)
                 contrasts(adf$OFscanner) <- 'contr.treatment'
                 # Inspect contrasts:
                 contrasts(adf$OFscanner)
                 
                 # ordered factoring of group
                 adf$OFgroup <- as.factor(adf$group)
                 # change factor to ordered factor:
                 adf$OFgroup <- as.ordered(adf$OFgroup)
                 # change contrast to treatment coding (difference curves)
                 contrasts(adf$OFgroup) <- 'contr.treatment'
                 # Inspect contrasts:
                 contrasts(adf$OFgroup)
                 
                 
                 # ordered factoring of med
                 adf$OFmed <- as.factor(adf$med)
                 # change factor to ordered factor:
                 adf$OFmed <- as.ordered(adf$OFmed)
                 # change contrast to treatment coding (difference curves)
                 contrasts(adf$OFmed) <- 'contr.treatment'
                 # Inspect contrasts:
                 contrasts(adf$OFmed)
                 
                 
                 
                 # Run the following gam models for each DV
              
                 # Null Model with covariates on DVs
                 gam1a <- gam(value  ~ 1 + SES_c + OFsex + meanFWD_c + etiv.long_c + OFscanner + OFmed + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth Age predictor
                 gam2a <- gam(value  ~ s(age_c,bs="cr",k=4) + SES_c + OFsex + meanFWD_c + etiv.long_c + OFscanner + OFmed + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 #Model with smooth age predictor and group
                 gam3a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFgroup + SES_c + OFsex + meanFWD_c + etiv.long_c + OFscanner + OFmed + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth age and group and interaction with age and group
                 gam4a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFgroup + s(age_c, by=OFgroup, bs='cr',k=4) + SES_c + OFsex + meanFWD_c + etiv.long_c + OFscanner + OFmed + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 
                 save(gam1a,gam2a,gam3a,gam4a, file = paste0("./FBA_gams_Total/LH/fullmodel/", dv_name,"_gams_models.Rdata"))
                 
                 # Model comparisons
                 
                 # 1 vs 2
                 # null model vs. smooth model
                 comp1 <- compareML(gam1a,gam2a)
                 table1 <- comp1$table
                 pvalue1 <- as.numeric(comp1$table[2,6])
                 table1$AIC <- comp1$AIC
                 table1$advice <- comp1$advice
                 
                 # 2 vs 3
                 # smooth model vs. main effect (sex) model
                 comp2 <- compareML(gam2a,gam3a)
                 table2 <- comp2$table
                 pvalue2 <- as.numeric(comp2$table[2,6])
                 table2$AIC <- comp2$AIC
                 table2$advice <- comp2$advice
                 
                 # 3 vs 4
                 #  main effect (sex) model vs. interaction model
                 comp3 <- compareML(gam3a,gam4a)
                 table3 <- comp3$table
                 pvalue3 <- as.numeric(comp3$table[2,6])
                 table3$AIC <- comp3$AIC
                 table3$advice <- comp3$advice

                 
                 #Saving output
                 
                 # save null model
                 ptable <- as.data.frame(summary(gam1a)$p.table) %>% rownames_to_column()
                 stable <- as.data.frame(summary(gam1a)$s.table) %>% rownames_to_column()
                 
                 # save smooth model
                 ptable2 <- as.data.frame(summary(gam2a)$p.table) %>% rownames_to_column()
                 stable2 <- as.data.frame(summary(gam2a)$s.table) %>% rownames_to_column()
                 
                 # save main effect (sex) model
                 ptable3 <- as.data.frame(summary(gam3a)$p.table) %>% rownames_to_column()
                 stable3 <- as.data.frame(summary(gam3a)$s.table) %>% rownames_to_column()
                 
                 # save interaction model
                 ptable4 <- as.data.frame(summary(gam4a)$p.table) %>% rownames_to_column()
                 stable4 <- as.data.frame(summary(gam4a)$s.table) %>% rownames_to_column()
                 
                 
                 
                 #Plotting GAMS

                 #plot trajectory
                 age_c <- round(seq(min(adf$age_c),max(adf$age_c),by=0.1),2)

                 lab <- gsub("value","",dv_name)
                 
                 if (pvalue1 < 0.05) {
                   if (pvalue2 < 0.05) {
                     if (pvalue3 < 0.05) {
                       
                       best_model = gam4a
                       
                       #predictor for gam4a (interaction)
                       plotdf3 <- data.frame(age_c=rep(age_c,2),OFgroup=c(rep("ADHD",length(age_c)),rep("CONTROL",length(age_c))),SID=c(rep("86",length(age_c)),rep("307",length(age_c))), SES_c = 0, OFsex = "female", meanFWD_c = 0, etiv.long_c = 0, OFmed = 0, OFscanner = "Tim_Trio")
                       pred3 <- predict(gam4a,plotdf3,exclude='s(SID)',se.fit=T)
                       plotdf3$pred <- pred3$fit
                       plotdf3$se <- pred3$se
                       plotdf3$lower <- plotdf3$pred - (1.96*(plotdf3$se))
                       plotdf3$upper <- plotdf3$pred + (1.96*(plotdf3$se))
                       plotdf3$age <- round(plotdf3$age_c + mean(adf$age),1)
                       plotdf3$OFgroup <- factor(plotdf3$OFgroup, levels=c("ADHD","CONTROL"))
                       
                       pred3a <- transform(cbind(data.frame(pred3),plotdf3))
                       colnames(pred3a)[1] <- "dv"
                       
                       
                       #Plot for gam4a (interaction)
                       best_plot <- ggplot() +
                         geom_smooth(data=pred3a,aes(age,dv,colour=OFgroup),method="gam",formula = y~s(x,k=4)) +
                         geom_ribbon(data=pred3a,aes(x=age,ymin=lower,ymax=upper,colour=OFgroup,fill=OFgroup),alpha=0.2)+
                         scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                         geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=2,alpha=0.3) +
                         geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=.3,alpha=0.3) +
                         scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                         theme_bw() +
                         theme_minimal(base_size = 12, base_family = "Arial") +
                         theme(axis.line = element_line(colour = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               legend.position="bottom") + 
                         ylab(lab) +
                         xlab("Age")
                       
                    }
                     else { best_model = gam3a
                     
                     #predictor for model gam3a (group)
                     plotdf2 <- data.frame(age_c=rep(age_c,2),OFgroup=c(rep("ADHD",length(age_c)),rep("CONTROL",length(age_c))),SID=c(rep("86",length(age_c)),rep("307",length(age_c))), SES_c = 0, OFsex = "female", meanFWD_c = 0, etiv.long_c = 0, Ofmed = 0, OFscanner = "Tim_Trio")
                     pred2 <- predict(gam3a,plotdf2,exclude='s(SID)',se.fit=T)
                     plotdf2$pred <- pred2$fit
                     plotdf2$se <- pred2$se
                     plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                     plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                     plotdf2$age <- round(plotdf2$age_c + mean(adf$age),1)
                     plotdf2$OFgroup <- factor(plotdf2$OFgroup, levels=c("ADHD","CONTROL"))
                     
                     pred2a <- transform(cbind(data.frame(pred2),plotdf2))
                     colnames(pred2a)[1] <- "dv"
                     
                     
                     #Plot for gam3a (group)
                     best_plot <- ggplot() +
                       geom_smooth(data=pred2a,aes(age,dv,colour=OFgroup),method="gam",formula = y~s(x,k=4)) +
                       geom_ribbon(data=pred2a,aes(x=age,ymin=lower,ymax=upper,colour=OFgroup,fill=OFgroup),alpha=0.2)+
                       scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                       geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=2,alpha=0.3) +
                       geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=.3,alpha=0.3) +
                       scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                       theme_bw() +
                       theme_minimal(base_size = 12, base_family = "Arial") +
                       theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             legend.position="bottom") + 
                       ylab(lab) +
                       facet_wrap(~OFgroup) +
                       xlab("Age") 
                     }
                     
                     } else { best_model = gam2a
                   
                   #predictor for model gam2a (Age)
                   plotdf <- data.frame(age_c=rep(age_c,2),OFgroup=c(rep("ADHD",length(age_c)),rep("CONTROL",length(age_c))),SID=c(rep("86",length(age_c)),rep("307",length(age_c))), SES_c = 0, OFsex = "female", meanFWD_c = 0, etiv.long_c = 0, OFmed = 0, OFscanner = "Tim_Trio")
                   pred1 <- predict(gam2a,plotdf,exclude='s(SID)',se.fit=T)
                   plotdf$pred <- pred1$fit
                   plotdf$se <- pred1$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$age <- round(plotdf$age_c + mean(adf$age),1)
                   plotdf$OFgroup <- factor(plotdf$OFgroup, levels=c("ADHD","CONTROL"))
                   
                   # averaging over group
                   pred_grp <- transform(cbind(data.frame(pred1),plotdf))
                   colnames(pred_grp)[1] <- "dv"
                   predicted_data_ADHD <- pred_grp[pred_grp$OFgroup=="ADHD",]
                   predicted_data_CONTROL <- pred_grp[pred_grp$OFgroup=="CONTROL",]
                   predicted_data_average <- predicted_data_ADHD
                   predicted_data_average$dv <- (predicted_data_ADHD$dv + predicted_data_CONTROL$dv)/2
                   predicted_data_average$upper <- (predicted_data_ADHD$upper + predicted_data_CONTROL$upper)/2
                   predicted_data_average$lower <- (predicted_data_ADHD$lower + predicted_data_CONTROL$lower)/2


                   #Plot for gam2a (Age)
                   best_plot <- ggplot() +
                     geom_smooth(data=predicted_data_average,aes(age,dv),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=predicted_data_average,aes(x=age,ymin=lower,ymax=upper),alpha=0.2)+
                     scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                     geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=2,alpha=0.3) +
                     scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("ADHD", "CONTROL")) +
                     geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFgroup),size=.3,alpha=0.3) +
                     theme_bw() +
                     theme_minimal(base_size = 12, base_family = "Arial") +
                     theme(axis.line = element_line(colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank(),
                           legend.position="bottom") + 
                     ylab(lab) +
                     xlab("Age")
                   }
                   
                   ggsave(filename=paste0("./FBA_gams_Total/LH/fullmodel/", dv_name, '.png'), plot=best_plot,width=6,height=5)
                   
                   tiff(paste0("./FBA_gams_Total/LH/fullmodel/", dv_name, "_check.tiff"), width = 9, height = 6, units = 'in', res = 150)
                   par(mfrow=c(2,2))
                   gam.check(best_model)
                   dev.off()
                   
                 }
                 
                 #combine all tables
                 table <- rbind.fill(table1,table2,table3,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                 table$dv <- dv_name
                 table
               })

LH_FBA_models_EXG3 <- rbindlist(models,fill=TRUE)
View(LH_FBA_models_EXG3)
write.csv(LH_FBA_models_EXG3,file= "./FBA_gams_Total/LH/fullmodel/lh_fba_models.csv",row.names=F)


LH_modelsFDR <- LH_FBA_models_EXG3[,c(1,6,7,10,14,18,19)]
LH_modelsFDR$p.fdr <- p.adjust(LH_modelsFDR$`Pr(>|t|)`,method="fdr")
LH_modelsFDR$pvalue.fdr <- p.adjust(LH_modelsFDR$`p.value`,method="fdr")
LH_modelsFDR$`p-value.fdr` <- p.adjust(LH_modelsFDR$`p-value`,method="fdr")
write_csv(LH_modelsFDR,path = "./FBA_gams_Total/LH/fullmodel/lh_fba_models_FDR.csv")

#######
