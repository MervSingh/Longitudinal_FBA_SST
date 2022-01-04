
# Study 2 CONTROLS - RIGHT HEMISPHERIC FRONTO-BASAL-GANGLIA FBA TRAJECTORIES

dir.create("./FBA_gams_Controls_only/RH")

# MODELS WITH COVARS (FULL MODEL)

# GAM models for EXG3 parameters
dir.create("./FBA_gams_Controls_only/RH/fullmodel")

models<-lapply(X=as.character(as.list(rh_dvs_CON)),
               df=rh_xlong_CON,
               FUN=function(dv_name, df) {
                 print(dv_name)
                 
                 # now create a new dataframe within the loop that is filtered to dv_name
                 adf<-df %>% filter(rh_dvs_CON==dv_name) %>%
                 mutate(SID = as.factor(SID),
                        Wave = as.factor(Wave), 
                        scanner = as.factor(scanner),
                        etiv.long_c = as.numeric(etiv.long_c),
                        meanFWD_c = as.numeric(meanFWD_c),
                        age_c = as.numeric(age_c),
                        sex = as.factor(sex),
                        hand = as.factor(hand), 
                        med = as.factor(med), 
                        SES_c = as.numeric(SES_c))
  
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
                 
                 
                 # Run the following gam models for each DV
              
                 # Null Model with covariates on DVs
                 gam1a <- gam(value  ~ 1 + SES_c + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth Age predictor
                 gam2a <- gam(value  ~ s(age_c,bs="cr",k=4) + SES_c + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 #Model with smooth age predictor and sex
                 gam3a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + SES_c + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth age and sex and interaction with age and group
                 gam4a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + s(age_c, by=OFsex, bs='cr',k=4) + SES_c + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 
                 save(gam1a,gam2a,gam3a,gam4a, file = paste0("./FBA_gams_Controls_only/RH/fullmodel/", dv_name,"_gams_models.Rdata"))
                 
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
                       plotdf3 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                       pred3 <- predict(gam4a,plotdf3,exclude='s(SID)',se.fit=T)
                       plotdf3$pred <- pred3$fit
                       plotdf3$se <- pred3$se
                       plotdf3$lower <- plotdf3$pred - (1.96*(plotdf3$se))
                       plotdf3$upper <- plotdf3$pred + (1.96*(plotdf3$se))
                       plotdf3$age <- round(plotdf3$age_c + mean(adf$age),1)
                       plotdf3$OFsex <- factor(plotdf3$OFsex, levels=c("female","male"))
                       
                       
                       pred3a <- transform(cbind(data.frame(pred3),plotdf3))
                       colnames(pred3a)[1] <- "dv"
                       
                       
                       #Plot for gam4a (interaction)
                       best_plot <- ggplot() +
                         geom_smooth(data=pred3a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                         geom_ribbon(data=pred3a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                         scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) + 
                         geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                         geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                         scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
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
                     
                     #predictor for model gam3a (sex)
                     plotdf2 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                     pred2 <- predict(gam3a,plotdf2,exclude='s(SID)',se.fit=T)
                     plotdf2$pred <- pred2$fit
                     plotdf2$se <- pred2$se
                     plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                     plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                     plotdf2$age <- round(plotdf2$age_c + mean(adf$age),1)
                     plotdf2$OFsex <- factor(plotdf2$OFsex, levels=c("female","male"))
                     
                     pred2a <- transform(cbind(data.frame(pred2),plotdf2))
                     colnames(pred2a)[1] <- "dv"
                     
                     
                     #Plot for gam3a (Sex)
                     best_plot <- ggplot() +
                       geom_smooth(data=pred2a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                       geom_ribbon(data=pred2a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                       scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                       geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                       scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       theme_bw() +
                       theme_minimal(base_size = 12, base_family = "Arial") +
                       theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             legend.position="bottom") + 
                       ylab(lab) +
                       facet_wrap(~OFsex) +
                       xlab("Age") 
                     }
                     
                     } else { best_model = gam2a
                   
                   #predictor for model gam2a (Age)
                   plotdf <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                   pred1 <- predict(gam2a,plotdf,exclude='s(SID)',se.fit=T)
                   plotdf$pred <- pred1$fit
                   plotdf$se <- pred1$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$age <- round(plotdf$age_c + mean(adf$age),1)
                   plotdf$OFsex <- factor(plotdf$OFsex, levels=c("female","male"))
                   
                   # #averaging over sex
                   pred_sex <- transform(cbind(data.frame(pred1),plotdf))
                   colnames(pred_sex)[1] <- "dv"
                   predicted_data_male <- pred_sex[pred_sex$OFsex=="male",]
                   predicted_data_female <- pred_sex[pred_sex$OFsex=="female",]
                   predicted_data_average <- predicted_data_male
                   predicted_data_average$dv <- (predicted_data_male$dv + predicted_data_female$dv)/2
                   predicted_data_average$upper <- (predicted_data_male$upper + predicted_data_female$upper)/2
                   predicted_data_average$lower <- (predicted_data_male$lower + predicted_data_female$lower)/2

                   
                   #Plot for gam2a (Age)
                   best_plot <- ggplot() +
                     geom_smooth(data=predicted_data_average,aes(age,dv),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=predicted_data_average,aes(x=age,ymin=lower,ymax=upper),alpha=0.2)+
                     scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                     scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
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
                   
                   ggsave(filename=paste0("./FBA_gams_Controls_only/RH/fullmodel/", dv_name, '.png'), plot=best_plot,width=6,height=5)
                   
                   tiff(paste0("./FBA_gams_Controls_only/RH/fullmodel/", dv_name, "_check.tiff"), width = 9, height = 6, units = 'in', res = 150)
                   par(mfrow=c(2,2))
                   gam.check(best_model)
                   dev.off()
                   
                 }
                 
                 #combine all tables
                 table <- rbind.fill(table1,table2,table3,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                 table$dv <- dv_name
                 table
               })

RH_FBA_models_EXG3 <- rbindlist(models,fill=TRUE)
View(RH_FBA_models_EXG3)
write.csv(RH_FBA_models_EXG3,file= "./FBA_gams_Controls_only/RH/fullmodel/rh_fba_models.csv",row.names=F)


RH_modelsFDR <- RH_FBA_models_EXG3[,c(1,6,7,10,14,18,19)]
RH_modelsFDR$p.fdr <- p.adjust(RH_modelsFDR$`Pr(>|t|)`,method="fdr")
RH_modelsFDR$pvalue.fdr <- p.adjust(RH_modelsFDR$`p.value`,method="fdr")
RH_modelsFDR$`p-value.fdr` <- p.adjust(RH_modelsFDR$`p-value`,method="fdr")
write_csv(RH_modelsFDR,path = "./FBA_gams_Controls_only/RH/fullmodel/rh_fba_models_FDR.csv")


###


# MODELS WITHOUT SES

# GAM models for EXG3 parameters
dir.create("./FBA_gams_Controls_only/RH/noSES")

models_noSES<-lapply(X=as.character(as.list(rh_dvs_CON)),
               df=rh_xlong_CON,
               FUN=function(dv_name, df) {
                 print(dv_name)

                 # now create a new dataframe within the loop that is filtered to dv_name
                 adf<-df %>% filter(rh_dvs_CON==dv_name) %>%
                   mutate(SID = as.factor(SID),
                          Wave = as.factor(Wave),
                          scanner = as.factor(scanner),
                          etiv.long_c = as.numeric(etiv.long_c),
                          meanFWD_c = as.numeric(meanFWD_c),
                          age_c = as.numeric(age_c),
                          sex = as.factor(sex),
                          hand = as.factor(hand),
                          med = as.factor(med),
                          SES_c = as.numeric(SES_c))

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


                 # Run the following gam models for each DV

                 # Null Model with covariates on DVs
                 gam1a <- gam(value  ~ 1 + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth Age predictor
                 gam2a <- gam(value  ~ s(age_c,bs="cr",k=4) + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 #Model with smooth age predictor and sex
                 gam3a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth age and sex and interaction with age and group
                 gam4a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + s(age_c, by=OFsex, bs='cr',k=4) + meanFWD_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option

                 save(gam1a,gam2a,gam3a,gam4a, file = paste0("./FBA_gams_Controls_only/RH/noSES/", dv_name,"_gams_models_noSES.Rdata"))

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
                       plotdf3 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                       pred3 <- predict(gam4a,plotdf3,exclude='s(SID)',se.fit=T)
                       plotdf3$pred <- pred3$fit
                       plotdf3$se <- pred3$se
                       plotdf3$lower <- plotdf3$pred - (1.96*(plotdf3$se))
                       plotdf3$upper <- plotdf3$pred + (1.96*(plotdf3$se))
                       plotdf3$age <- round(plotdf3$age_c + mean(adf$age),1)
                       plotdf3$OFsex <- factor(plotdf3$OFsex, levels=c("female","male"))
                       
                       pred3a <- transform(cbind(data.frame(pred3),plotdf3))
                       colnames(pred3a)[1] <- "dv"
                       

                       #Plot for gam4a (interaction)
                       best_plot <- ggplot() +
                         geom_smooth(data=pred3a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                         geom_ribbon(data=pred3a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                         scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                         geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                         geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                         scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
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

                     #predictor for model gam3a (sex)
                     plotdf2 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                     pred2 <- predict(gam3a,plotdf2,exclude='s(SID)',se.fit=T)
                     plotdf2$pred <- pred2$fit
                     plotdf2$se <- pred2$se
                     plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                     plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                     plotdf2$age <- round(plotdf2$age_c + mean(adf$age),1)
                     plotdf2$OFsex <- factor(plotdf2$OFsex, levels=c("female","male"))
                     
                     pred2a <- transform(cbind(data.frame(pred2),plotdf2))
                     colnames(pred2a)[1] <- "dv"
                     

                     #Plot for gam3a (Sex)
                     best_plot <- ggplot() +
                       geom_smooth(data=pred2a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                       geom_ribbon(data=pred2a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                       scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                       geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                       scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       theme_bw() +
                       theme_minimal(base_size = 12, base_family = "Arial") +
                       theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             legend.position="bottom") +
                       ylab(lab) +
                       facet_wrap(~OFsex) +
                       xlab("Age")
                     }

                   } else { best_model = gam2a

                   #predictor for model gam2a (Age)
                   plotdf <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), meanFWD_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                   pred1 <- predict(gam2a,plotdf,exclude='s(SID)',se.fit=T)
                   plotdf$pred <- pred1$fit
                   plotdf$se <- pred1$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$age <- round(plotdf$age_c + mean(adf$age),1)
                   plotdf$OFsex <- factor(plotdf$OFsex, levels=c("female","male"))

                   # #averaging over sex
                   pred_sex <- transform(cbind(data.frame(pred1),plotdf))
                   colnames(pred_sex)[1] <- "dv"
                   predicted_data_male <- pred_sex[pred_sex$OFsex=="male",]
                   predicted_data_female <- pred_sex[pred_sex$OFsex=="female",]
                   predicted_data_average <- predicted_data_male
                   predicted_data_average$dv <- (predicted_data_male$dv + predicted_data_female$dv)/2
                   predicted_data_average$upper <- (predicted_data_male$upper + predicted_data_female$upper)/2
                   predicted_data_average$lower <- (predicted_data_male$lower + predicted_data_female$lower)/2


                   #Plot for gam2a (Age)
                   best_plot <- ggplot() +
                     geom_smooth(data=predicted_data_average,aes(age,dv),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=predicted_data_average,aes(x=age,ymin=lower,ymax=upper),alpha=0.2)+
                     scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                     scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
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

                   ggsave(filename=paste0("./FBA_gams_Controls_only/RH/noSES/", dv_name, '_noSES.png'), plot=best_plot,width=6,height=5)

                   tiff(paste0("./FBA_gams_Controls_only/RH/noSES/", dv_name, "_check_noSES.tiff"), width = 9, height = 6, units = 'in', res = 150)
                   par(mfrow=c(2,2))
                   gam.check(best_model)
                   dev.off()

                 }

                 #combine all tables
                 table <- rbind.fill(table1,table2,table3,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                 table$dv <- dv_name
                 table
               })

RH_FBA_models_EXG3_noSES <- rbindlist(models_noSES,fill=TRUE)
View(RH_FBA_models_EXG3_noSES)
write.csv(RH_FBA_models_EXG3_noSES,file= "./FBA_gams_Controls_only/RH/noSES/rh_fba_models_noSES.csv",row.names=F)


RH_modelsFDR_noSES <- RH_FBA_models_EXG3_noSES[,c(1,6,7,10,14,18,19)]
RH_modelsFDR_noSES$p.fdr <- p.adjust(RH_modelsFDR_noSES$`Pr(>|t|)`,method="fdr")
RH_modelsFDR_noSES$pvalue.fdr <- p.adjust(RH_modelsFDR_noSES$`p.value`,method="fdr")
RH_modelsFDR_noSES$`p-value.fdr` <- p.adjust(RH_modelsFDR_noSES$`p-value`,method="fdr")
write_csv(RH_modelsFDR_noSES,path = "./FBA_gams_Controls_only/RH/noSES/rh_fba_models_FDR_noSES.csv")



###



# MODELS WITHOUT FWD

# GAM models for EXG3 parameters
dir.create("./FBA_gams_Controls_only/RH/noFWD")

models_noFWD<-lapply(X=as.character(as.list(rh_dvs_CON)),
               df=rh_xlong_CON,
               FUN=function(dv_name, df) {
                 print(dv_name)

                 # now create a new dataframe within the loop that is filtered to dv_name
                 adf<-df %>% filter(rh_dvs_CON==dv_name) %>%
                   mutate(SID = as.factor(SID),
                          Wave = as.factor(Wave),
                          scanner = as.factor(scanner),
                          etiv.long_c = as.numeric(etiv.long_c),
                          meanFWD_c = as.numeric(meanFWD_c),
                          age_c = as.numeric(age_c),
                          sex = as.factor(sex),
                          hand = as.factor(hand),
                          med = as.factor(med),
                          SES_c = as.numeric(SES_c))

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


                 # Run the following gam models for each DV

                 # Null Model with covariates on DVs
                 gam1a <- gam(value  ~ 1 + SES_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth Age predictor
                 gam2a <- gam(value  ~ s(age_c,bs="cr",k=4) + SES_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 #Model with smooth age predictor and sex
                 gam3a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + SES_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth age and sex and interaction with age and group
                 gam4a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + s(age_c, by=OFsex, bs='cr',k=4) + SES_c + etiv.long_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option

                 save(gam1a,gam2a,gam3a,gam4a, file = paste0("./FBA_gams_Controls_only/RH/noFWD/", dv_name,"_gams_models_noFWD.Rdata"))

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
                       plotdf3 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                       pred3 <- predict(gam4a,plotdf3,exclude='s(SID)',se.fit=T)
                       plotdf3$pred <- pred3$fit
                       plotdf3$se <- pred3$se
                       plotdf3$lower <- plotdf3$pred - (1.96*(plotdf3$se))
                       plotdf3$upper <- plotdf3$pred + (1.96*(plotdf3$se))
                       plotdf3$age <- round(plotdf3$age_c + mean(adf$age),1)
                       plotdf3$OFsex <- factor(plotdf3$OFsex, levels=c("female","male"))
                       
                       pred3a <- transform(cbind(data.frame(pred3),plotdf3))
                       colnames(pred3a)[1] <- "dv"
                       

                       #Plot for gam4a (interaction)
                       best_plot <- ggplot() +
                         geom_smooth(data=pred3a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                         geom_ribbon(data=pred3a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                         scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                         geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                         geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                         scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
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

                     #predictor for model gam3a (sex)
                     plotdf2 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                     pred2 <- predict(gam3a,plotdf2,exclude='s(SID)',se.fit=T)
                     plotdf2$pred <- pred2$fit
                     plotdf2$se <- pred2$se
                     plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                     plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                     plotdf2$age <- round(plotdf2$age_c + mean(adf$age),1)
                     plotdf2$OFsex <- factor(plotdf2$OFsex, levels=c("female","male"))
                     
                     pred2a <- transform(cbind(data.frame(pred2),plotdf2))
                     colnames(pred2a)[1] <- "dv"
                     

                     #Plot for gam3a (Sex)
                     best_plot <- ggplot() +
                       geom_smooth(data=pred2a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                       geom_ribbon(data=pred2a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                       scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                       geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                       scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       theme_bw() +
                       theme_minimal(base_size = 12, base_family = "Arial") +
                       theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             legend.position="bottom") +
                       ylab(lab) +
                       facet_wrap(~OFsex) +
                       xlab("Age")
                     }

                   } else { best_model = gam2a

                   #predictor for model gam2a (Age)
                   plotdf <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, etiv.long_c = 0, OFscanner = "Tim_Trio")
                   pred1 <- predict(gam2a,plotdf,exclude='s(SID)',se.fit=T)
                   plotdf$pred <- pred1$fit
                   plotdf$se <- pred1$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$age <- round(plotdf$age_c + mean(adf$age),1)
                   plotdf$OFsex <- factor(plotdf$OFsex, levels=c("female","male"))

                   # #averaging over sex
                   pred_sex <- transform(cbind(data.frame(pred1),plotdf))
                   colnames(pred_sex)[1] <- "dv"
                   predicted_data_male <- pred_sex[pred_sex$OFsex=="male",]
                   predicted_data_female <- pred_sex[pred_sex$OFsex=="female",]
                   predicted_data_average <- predicted_data_male
                   predicted_data_average$dv <- (predicted_data_male$dv + predicted_data_female$dv)/2
                   predicted_data_average$upper <- (predicted_data_male$upper + predicted_data_female$upper)/2
                   predicted_data_average$lower <- (predicted_data_male$lower + predicted_data_female$lower)/2


                   #Plot for gam2a (Age)
                   best_plot <- ggplot() +
                     geom_smooth(data=predicted_data_average,aes(age,dv),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=predicted_data_average,aes(x=age,ymin=lower,ymax=upper),alpha=0.2)+
                     scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                     scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
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

                   ggsave(filename=paste0("./FBA_gams_Controls_only/RH/noFWD/", dv_name, '_noFWD.png'), plot=best_plot,width=6,height=5)

                   tiff(paste0("./FBA_gams_Controls_only/RH/noFWD/", dv_name, "_check_noFWD.tiff"), width = 9, height = 6, units = 'in', res = 150)
                   par(mfrow=c(2,2))
                   gam.check(best_model)
                   dev.off()

                 }

                 #combine all tables
                 table <- rbind.fill(table1,table2,table3,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                 table$dv <- dv_name
                 table
               })

RH_FBA_models_EXG3_noFWD <- rbindlist(models_noFWD,fill=TRUE)
View(RH_FBA_models_EXG3_noFWD)
write.csv(RH_FBA_models_EXG3_noFWD,file= "./FBA_gams_Controls_only/RH/noFWD/rh_fba_models_noFWD.csv",row.names=F)


RH_modelsFDR_noFWD <- RH_FBA_models_EXG3_noFWD[,c(1,6,7,10,14,18,19)]
RH_modelsFDR_noFWD$p.fdr <- p.adjust(RH_modelsFDR_noFWD$`Pr(>|t|)`,method="fdr")
RH_modelsFDR_noFWD$pvalue.fdr <- p.adjust(RH_modelsFDR_noFWD$`p.value`,method="fdr")
RH_modelsFDR_noFWD$`p-value.fdr` <- p.adjust(RH_modelsFDR_noFWD$`p-value`,method="fdr")
write_csv(RH_modelsFDR_noFWD,path = "./FBA_gams_Controls_only/RH/noFWD/rh_fba_models_FDR_noFWD.csv")



###



# MODELS WITHOUT ICV

# GAM models for EXG3 parameters
dir.create("./FBA_gams_Controls_only/RH/noICV")

models_noICV<-lapply(X=as.character(as.list(rh_dvs_CON)),
               df=rh_xlong_CON,
               FUN=function(dv_name, df) {
                 print(dv_name)

                 # now create a new dataframe within the loop that is filtered to dv_name
                 adf<-df %>% filter(rh_dvs_CON==dv_name) %>%
                   mutate(SID = as.factor(SID),
                          Wave = as.factor(Wave),
                          scanner = as.factor(scanner),
                          etiv.long_c = as.numeric(etiv.long_c),
                          meanFWD_c = as.numeric(meanFWD_c),
                          age_c = as.numeric(age_c),
                          sex = as.factor(sex),
                          hand = as.factor(hand),
                          med = as.factor(med),
                          SES_c = as.numeric(SES_c))

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


                 # Run the following gam models for each DV

                 # Null Model with covariates on DVs
                 gam1a <- gam(value  ~ 1 + SES_c + meanFWD_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth Age predictor
                 gam2a <- gam(value  ~ s(age_c,bs="cr",k=4) + SES_c + meanFWD_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 #Model with smooth age predictor and sex
                 gam3a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + SES_c + meanFWD_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option
                 # Model with smooth age and sex and interaction with age and group
                 gam4a <- gam(value  ~ s(age_c,bs="cr",k=4) + OFsex + s(age_c, by=OFsex, bs='cr',k=4) + SES_c + meanFWD_c + OFscanner + s(SID,bs='re'), data=adf, method = "ML") # add 'family' option

                 save(gam1a,gam2a,gam3a,gam4a, file = paste0("./FBA_gams_Controls_only/RH/noICV/", dv_name,"_gams_models_noICV.Rdata"))

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
                       plotdf3 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, OFscanner = "Tim_Trio")
                       pred3 <- predict(gam4a,plotdf3,exclude='s(SID)',se.fit=T)
                       plotdf3$pred <- pred3$fit
                       plotdf3$se <- pred3$se
                       plotdf3$lower <- plotdf3$pred - (1.96*(plotdf3$se))
                       plotdf3$upper <- plotdf3$pred + (1.96*(plotdf3$se))
                       plotdf3$age <- round(plotdf3$age_c + mean(adf$age),1)
                       plotdf3$OFsex <- factor(plotdf3$OFsex, levels=c("female","male"))
                       
                       pred3a <- transform(cbind(data.frame(pred3),plotdf3))
                       colnames(pred3a)[1] <- "dv"
                       

                       #Plot for gam4a (interaction)
                       best_plot <- ggplot() +
                         geom_smooth(data=pred3a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                         geom_ribbon(data=pred3a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                         scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                         geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                         geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                         scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
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

                     #predictor for model gam3a (sex)
                     plotdf2 <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, OFscanner = "Tim_Trio")
                     pred2 <- predict(gam3a,plotdf2,exclude='s(SID)',se.fit=T)
                     plotdf2$pred <- pred2$fit
                     plotdf2$se <- pred2$se
                     plotdf2$lower <- plotdf2$pred - (1.96*(plotdf2$se))
                     plotdf2$upper <- plotdf2$pred + (1.96*(plotdf2$se))
                     plotdf2$age <- round(plotdf2$age_c + mean(adf$age),1)
                     plotdf2$OFsex <- factor(plotdf2$OFsex, levels=c("female","male"))
                     
                     pred2a <- transform(cbind(data.frame(pred2),plotdf2))
                     colnames(pred2a)[1] <- "dv"
                     

                     #Plot for gam3a (Sex)
                     best_plot <- ggplot() +
                       geom_smooth(data=pred2a,aes(age,dv,colour=OFsex),method="gam",formula = y~s(x,k=4)) +
                       geom_ribbon(data=pred2a,aes(x=age,ymin=lower,ymax=upper,colour=OFsex,fill=OFsex),alpha=0.2)+
                       scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                       geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
                       scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                       theme_bw() +
                       theme_minimal(base_size = 12, base_family = "Arial") +
                       theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(),
                             legend.position="bottom") +
                       ylab(lab) +
                       facet_wrap(~OFsex) +
                       xlab("Age")
                     }

                   } else { best_model = gam2a

                   #predictor for model gam2a (Age)
                   plotdf <- data.frame(age_c=rep(age_c,2),OFsex=c(rep("female",length(age_c)),rep("male",length(age_c))),SID=c(rep("85",length(age_c)),rep("15",length(age_c))), SES_c = 0, meanFWD_c = 0, OFscanner = "Tim_Trio")
                   pred1 <- predict(gam2a,plotdf,exclude='s(SID)',se.fit=T)
                   plotdf$pred <- pred1$fit
                   plotdf$se <- pred1$se
                   plotdf$lower <- plotdf$pred - (1.96*(plotdf$se))
                   plotdf$upper <- plotdf$pred + (1.96*(plotdf$se))
                   plotdf$age <- round(plotdf$age_c + mean(adf$age),1)
                   plotdf$OFsex <- factor(plotdf$OFsex, levels=c("female","male"))

                   # #averaging over sex
                   pred_sex <- transform(cbind(data.frame(pred1),plotdf))
                   colnames(pred_sex)[1] <- "dv"
                   predicted_data_male <- pred_sex[pred_sex$OFsex=="male",]
                   predicted_data_female <- pred_sex[pred_sex$OFsex=="female",]
                   predicted_data_average <- predicted_data_male
                   predicted_data_average$dv <- (predicted_data_male$dv + predicted_data_female$dv)/2
                   predicted_data_average$upper <- (predicted_data_male$upper + predicted_data_female$upper)/2
                   predicted_data_average$lower <- (predicted_data_male$lower + predicted_data_female$lower)/2


                   #Plot for gam2a (Age)
                   best_plot <- ggplot() +
                     geom_smooth(data=predicted_data_average,aes(age,dv),method="gam",formula = y~s(x,k=4)) +
                     geom_ribbon(data=predicted_data_average,aes(x=age,ymin=lower,ymax=upper),alpha=0.2)+
                     scale_fill_manual(name = "Fitted model", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_point(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=2,alpha=0.3) +
                     scale_colour_manual(name = "Raw data", values=c("#F8766D","#00BFC4"), labels = c("male", "female")) +
                     geom_line(data=adf,aes(x=age,y=value,group=SID,colour=OFsex),size=.3,alpha=0.3) +
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

                   ggsave(filename=paste0("./FBA_gams_Controls_only/RH/noICV/", dv_name, '_noICV.png'), plot=best_plot,width=6,height=5)

                   tiff(paste0("./FBA_gams_Controls_only/RH/noICV/", dv_name, "_check_noICV.tiff"), width = 9, height = 6, units = 'in', res = 150)
                   par(mfrow=c(2,2))
                   gam.check(best_model)
                   dev.off()

                 }

                 #combine all tables
                 table <- rbind.fill(table1,table2,table3,ptable,stable,ptable2,stable2,ptable3,stable3,ptable4,stable4)
                 table$dv <- dv_name
                 table
               })

RH_FBA_models_EXG3_noICV <- rbindlist(models_noICV,fill=TRUE)
View(RH_FBA_models_EXG3_noICV)
write.csv(RH_FBA_models_EXG3_noICV,file= "./FBA_gams_Controls_only/RH/noICV/rh_fba_models_noICV.csv",row.names=F)


RH_modelsFDR_noICV <- RH_FBA_models_EXG3_noICV[,c(1,6,7,10,14,18,19)]
RH_modelsFDR_noICV$p.fdr <- p.adjust(RH_modelsFDR_noICV$`Pr(>|t|)`,method="fdr")
RH_modelsFDR_noICV$pvalue.fdr <- p.adjust(RH_modelsFDR_noICV$`p.value`,method="fdr")
RH_modelsFDR_noICV$`p-value.fdr` <- p.adjust(RH_modelsFDR_noICV$`p-value`,method="fdr")
write_csv(RH_modelsFDR_noICV,path = "./FBA_gams_Controls_only/RH/noICV/rh_fba_models_FDR_noICV.csv")



###
