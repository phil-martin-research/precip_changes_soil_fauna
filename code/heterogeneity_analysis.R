# script to explore heterogeneity in the effects of precipitation changes on soil and litter fauna 

rm(list = ls())

pacman::p_load(tidyverse,metafor,cowplot,orchaRd,ggbeeswarm,tidyr,ggthemes,sp,broom,lemon,MuMIn,glmulti,PerformanceAnalytics,GGally,gt,geodata,
               ggmap,mapproj,glmulti)

#read in .csv files with soil fauna data
abundance<- read_csv("data/abundance_data.csv")
richness<- read_csv("data/richness_data.csv")
shannon<- read_csv("data/shannon_data.csv")

###############################################################################
#1 - data tidying##############################################################
###############################################################################

# create a unit-level random effect to model residual variance in metafor
abundance$obsID <- 1:nrow(abundance)
# mean-centering year of publication to help with interpretation
abundance$year.c <- as.vector(scale(abundance$study_year, scale = F))

#complete cases for the variable about percentage change in precipitation and body size
abundance_complete<-abundance[complete.cases(abundance$perc_annual_dist,abundance$Functional_group_size.y,abundance$above_below,
                                             abundance$exoskeleton,abundance$obsID,abundance$year.c),]

# create a unit-level random effect to model residual variance in metafor
richness$obsID <- 1:nrow(richness)
# mean-centering year of publication to help with interpretation
richness$year.c <- as.vector(scale(richness$study_year, scale = F))

#we lose 21 comparisons for abundance
richness_complete<-richness[complete.cases(richness$perc_annual_dist,
                                           richness$Functional_group_size.y,
                                           richness$above_below,
                                           richness$exoskeleton,
                                           richness$obsID,
                                           richness$year.c),]
#we lose 5 comparisons for richness

# create a unit-level random effect to model residual variance in metafor
shannon$obsID <- 1:nrow(shannon)
# mean-centering year of publication to help with interpretation
shannon$year.c <- as.vector(scale(shannon$study_year, scale = F))

shannon_complete<-shannon[complete.cases(shannon$perc_annual_dist,
                                         shannon$Functional_group_size.y,
                                         shannon$above_below,
                                         shannon$exoskeleton,
                                         shannon$obsID,
                                         shannon$year.c),]
#we lose 2 comparisons for shannon wiener

#allow model evalution with MuMIn package
eval(metafor:::.MuMIn)

#incorporate data about location aridity
abundance_complete$arid_class<-ifelse(abundance_complete$aridity>0.65,"Humid","Dry")

###############################################################################
#1 - models of heterogeneity in response of abundance of soil and litter fauna#
###############################################################################


#test all models against each other

#1 - null model
M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#2 - impact of precipitation changes
M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist ,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#3 - impact of precipitation changes plus year effect
M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#4 - impact of precipitation changes plus small study effect
M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#5 - impact of precipitation changes plus year and small study effect
M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#6 - impacts of precipitation change are modified by whether organisms are above or belowground
M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#7 - impacts of precipitation change are modified by whether organisms are above or belowground plus year effect
M6<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#8 - impacts of precipitation change are modified by whether organisms are above or belowground plus small study effect
M7<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#9 - impacts of precipitation change are modified by whether organisms are above or belowground plus year and small study effect
M8<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#10 - precipitation change is modified by functional group size
M9<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#11 - precipitation change is modified by functional group size plus year effect
M10<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#12 - precipitation change is modified by functional group size plus small study effect
M11<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#13 - precipitation change is modified by functional group size plus year and small study
M12<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#14 - impacts of precipitation change is modified by whether organism has an exoskeleton or not
M13<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#15 - impacts of precipitation change is modified by whether organism has an exoskeleton or not plus year effect
M14<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#16 - impacts of precipitation change is modified by whether organism has an exoskeleton or not plus small study effect
M15<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#17 - impacts of precipitation change is modified by whether organism has an exoskeleton or not plus year and small study
M16<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#18 - impacts of precipitation change is modified by whether ecosystem is humid or dry
M17<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*arid_class+I(perc_annual_dist^2)*arid_class-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#19 - impacts of precipitation change is modified by whether ecosystem is humid or dry plus year effect
M18<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*arid_class+I(perc_annual_dist^2)*arid_class+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#20 - impacts of precipitation change is modified by whether ecosystem is humid or dry plus small study effect
M19<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*arid_class+I(perc_annual_dist^2)*arid_class+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)
#21 - impacts of precipitation change is modified by whether ecosystem is humid or dry plus year and small study
M20<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*arid_class+I(perc_annual_dist^2)*arid_class+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)

#create a model selection table of these models
abun_model_sel<-model.sel(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,M16,M17,M18,M19,M20)

#M12 is the best model that includes an interaction between functional group size and change in precipitation
#as well as correction for year effect and small study effect

#create dataframe for all models
abundance_model_sel_df<-data.frame(abun_model_sel)
abundance_model_sel_df$model<-row.names(abundance_model_sel_df)

#save this model table
abundance_model_sel_table<-abundance_model_sel_df%>%
  select(mods,df,logLik,AICc,delta,model)%>%
  mutate(across(where(is.numeric), round, 2))%>%
  gt()
#export this to a word file
abundance_model_sel_table%>%gtsave("figures/for_paper/abundance_selection_table.docx")

#subset to give top models
abun_model_sel_sub<-subset(abun_model_sel, delta <= 2, recalc.weights=FALSE)

#calculate model averaged coefficients
abun_coefs<-model.avg(abun_model_sel_sub)
summary_abun_coefs<-summary(abun_coefs)

#save the full average coefficient matrix
abun_coef_full<-data.frame(summary_abun_coefs$coefmat.full)

#tidy up this table
abun_coef_full$Variable<-row.names(abun_coef_full)
abun_coef_full_table<-abun_coef_full%>%
  relocate(Variable)%>%
  mutate(Estimate=round(Estimate,4),
         Std..Error=round(Std..Error,4),
         z.value=round(z.value,4),
         Pr...z..=round(Pr...z..,4))%>%
  rename(SE=Std..Error,
         p_value=Pr...z..)%>%
  remove_rownames()%>%
  mutate(Variable=str_replace(Variable,"Functional_group_size.y",""),
         Variable=str_replace(Variable,"perc_annual_dist","Precipitation change"),
         Variable=str_replace(Variable,"year.c","Decline effect"),
         Variable=str_replace(Variable,"sqrt_inv_n_tilda","Small study size"))%>%
  gt()
  
abun_coef_full_table%>%gtsave("figures/for_paper/abundance_coef_table.docx")


#copy and save the model formula of the best model
M12_formula<-(~perc_annual_dist*Functional_group_size.y+year.c+sqrt_inv_n_tilda-1)

#create dataframe with new data for predictions
new_data<-data.frame(expand.grid(
  perc_annual_dist = seq(min(abundance_complete$perc_annual_dist),max(abundance_complete$perc_annual_dist),0.1),
  Functional_group_size.y=levels(as.factor(abundance_complete$Functional_group_size.y)),
  year.c=mean(abundance_complete$year.c),
  sqrt_inv_n_tilda=mean(abundance_complete$sqrt_inv_n_tilda)))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M12_formula,data=new_data)

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(M12,newmods=predgrid))

#attach predictions to variables for plotting
new_data <- cbind(new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot data for model
#first subset to remove very large effect sizes
abundance_sub<-abundance_complete%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size.y,"microfauna","mesofauna","macrofauna"))

#remove predictions for microfauna and macrofauna where the changes in precipitation are outside the observed values
abundance_sub%>%
  group_by(Functional_group_size.y)%>%
  summarise(min_precip=min(perc_annual_dist),
            max_precip=max(perc_annual_dist))

new_data_micro<-subset(new_data,Functional_group_size.y=="microfauna"&perc_annual_dist<=100)
new_data_meso<-subset(new_data,Functional_group_size.y!="microfauna")

#merge these together
new_data_merge<-rbind(new_data_micro,new_data_meso)

#add data on sample size for each group
sample_size_label<-abundance_sub%>%
  group_by(Functional_group_size)%>%
  summarise(k=length(yi),study_n=n_distinct(Study_ID))%>%
  mutate(k_label=paste("k = ",k," (",study_n,")",sep = ""))

#new version of size and annual change plot
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size.y,"microfauna","mesofauna","macrofauna"))%>%
ggplot(aes(x=perc_annual_dist,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_rep_wrap(~Functional_group_size,repeat.tick.labels = TRUE)+
  geom_point(data=abundance_sub,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  xlab("Change in annual precipitation (%)")+
  ylab("Change in soil fauna abundance\n(log response ratio)")+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_fill_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_size_continuous(range = c(1,10),transform = "sqrt")+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)+
  geom_text(data=sample_size_label,aes(x=150,y=4,label=k_label),colour="black")

#save plot
ggsave("figures/Figure_3.pdf",width = 20,height = 8,units = "cm",dpi = 300)


#make predictions for decline effect
new_data_decline_ab<-data.frame(expand.grid(
  perc_annual_dist = mean(abundance_complete$perc_annual_dist),
  Functional_group_size.y=levels(as.factor(abundance_complete$Functional_group_size.y)),
  year.c=seq(min(abundance_complete$year.c),max(abundance_complete$year.c)),
  sqrt_inv_n_tilda=mean(abundance_complete$sqrt_inv_n_tilda)))

#create a model matrix and remove the intercept
predgrid_ab_decline<-model.matrix(M12_formula,data=new_data_decline_ab)

#predict onto the new model matrix
mypreds_ab_decline<-data.frame(predict.rma(M12,newmods=predgrid_ab_decline))

#attach predictions to variables for plotting
new_data_decline_ab <- cbind(new_data_decline_ab, mypreds_ab_decline[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot predictions for decline effect
new_data_decline_ab%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size.y,"microfauna","mesofauna","macrofauna"))%>%
  ggplot(aes(x=year.c+mean(abundance_sub$study_year),y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_rep_wrap(~Functional_group_size,repeat.tick.labels = TRUE)+
  geom_point(data=abundance_sub,aes(x=study_year,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  xlab("Year of publicaton")+
  ylab("Change in soil fauna abundance\n(log response ratio)")+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_fill_manual(values = c("#02475f","#c3386b","#e0b500"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)

#save plot
ggsave("figures/abundance_decline_effect.png",width = 20,height = 8,units = "cm",dpi = 300)

#make predictions for small study effect
new_data_study_size_ab<-data.frame(expand.grid(
  perc_annual_dist = mean(abundance_complete$perc_annual_dist),
  Functional_group_size.y=levels(as.factor(abundance_complete$Functional_group_size.y)),
  year.c=mean(abundance_complete$year.c),
  sqrt_inv_n_tilda=seq(min(abundance_complete$sqrt_inv_n_tilda),max(abundance_complete$sqrt_inv_n_tilda),0.01)))

#create a model matrix and remove the intercept
predgrid_ab_study<-model.matrix(M12_formula,data=new_data_study_size_ab)

#predict onto the new model matrix
mypreds_ab_study<-data.frame(predict.rma(M12,newmods=predgrid_ab_study))

#attach predictions to variables for plotting
new_data_study_size_ab <- cbind(new_data_study_size_ab, mypreds_ab_study[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot predictions for decline effect
new_data_study_size_ab%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size.y,"microfauna","mesofauna","macrofauna"))%>%
  ggplot(aes(x=sqrt_inv_n_tilda,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_rep_wrap(~Functional_group_size,repeat.tick.labels = TRUE)+
  geom_point(data=abundance_sub,aes(x=sqrt_inv_n_tilda,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  xlab("Study size")+
  ylab("Change in soil fauna abundance\n(log response ratio)")+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_fill_manual(values = c("#02475f","#c3386b","#e0b500"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)

#save plot
ggsave("figures/abundance_study_size_effect.png",width = 20,height = 8,units = "cm",dpi = 300)


#1.2 abundance models of different groups of mesofauna

#filter data and reformat to allow comparison between acari and collembola
mesofauna_abundance<-abundance_complete%>%
  filter(Functional_group_size.y=="mesofauna")%>%
  mutate(acari_collembola=if_else(Highest_taxonomic_resolution=="Mesostigmata"|
                                 Highest_taxonomic_resolution=="Oribatida"|
                                 Highest_taxonomic_resolution=="Acari"|
                                 Highest_taxonomic_resolution=="Prostigmata","Acari",
                         if_else(Highest_taxonomic_resolution=="Entomobryomorpha"|
                                 Highest_taxonomic_resolution=="Poduromorpha"|
                                 Highest_taxonomic_resolution=="Symphypleona"|
                                 Highest_taxonomic_resolution=="Neelipleona","Collembola",Highest_taxonomic_resolution)))%>%
  filter(acari_collembola=="Acari"|acari_collembola=="Collembola")



mesofauna_M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID/obsID,data=mesofauna_abundance)
mesofauna_M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*acari_collembola-1,random=~1|Study_ID/Site_ID/obsID,data=mesofauna_abundance)
mesofauna_M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*acari_collembola+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=mesofauna_abundance)
mesofauna_M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*acari_collembola+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=mesofauna_abundance)
mesofauna_M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*acari_collembola+year.c+sqrt_inv_n_tilda-1-1,random=~1|Study_ID/Site_ID/obsID,data=mesofauna_abundance)



#create a model selection table of these models
meso_abun_model_sel<-model.sel(mesofauna_M1,mesofauna_M2,mesofauna_M3,mesofauna_M4,mesofauna_M5)

#create dataframe for all models
meso_abundance_model_sel_df<-data.frame(meso_abun_model_sel)
meso_abundance_model_sel_df$model<-row.names(meso_abundance_model_sel_df)

#save this model table
meso_abundance_model_sel_table<-meso_abundance_model_sel_df%>%
  select(mods,df,logLik,AICc,delta,model)%>%
  mutate(across(where(is.numeric), round, 2))%>%
  gt()
#export this to a word file
meso_abundance_model_sel_table%>%gtsave("figures/meso_abundance_selection_table.docx")

#subset to give top models
meso_abun_model_sel_sub<-subset(meso_abun_model_sel, delta <= 2, recalc.weights=FALSE)

#calculate model averaged coefficients
meso_abun_coefs<-model.avg(meso_abun_model_sel_sub)
summary_meso_abun_coefs<-summary(meso_abun_coefs)

#save the full average coefficient matrix
meso_abun_coef_full<-data.frame(summary_meso_abun_coefs$coefmat.full)

#tidy up this table
meso_abun_coef_full$Variable<-row.names(meso_abun_coef_full)
meso_abun_coef_full_table<-meso_abun_coef_full%>%
  relocate(Variable)%>%
  mutate(Estimate=round(Estimate,4),
         Std..Error=round(Std..Error,4),
         z.value=round(z.value,4),
         Pr...z..=round(Pr...z..,4))%>%
  rename(SE=Std..Error,
         p_value=Pr...z..)%>%
  remove_rownames()%>%
  mutate(Variable=str_replace(Variable,"Functional_group_size.y",""),
         Variable=str_replace(Variable,"perc_annual_dist","Precipitation change"),
         Variable=str_replace(Variable,"year.c","Decline effect"),
         Variable=str_replace(Variable,"sqrt_inv_n_tilda","Small study size"),
         Variable=str_replace(Variable,"acari_collembolaAcari","Acari"),
         Variable=str_replace(Variable,"acari_collembolaCollembola","Collembola"),
         Variable=str_replace(Variable,"intrcpt","Intercept"))%>%
  gt()

meso_abun_coef_full_table%>%gtsave("figures/meso_abundance_coef_table.docx")

#produce a figure for the most parsimonious model

#copy and save the model formula of the best model
M4_formula<-(~perc_annual_dist*acari_collembola+sqrt_inv_n_tilda-1)

#create dataframe with new data for predictions
meso_new_data<-data.frame(expand.grid(
  perc_annual_dist = seq(min(mesofauna_abundance$perc_annual_dist),max(mesofauna_abundance$perc_annual_dist),0.1),
  acari_collembola=levels(as.factor(mesofauna_abundance$acari_collembola)),
  year.c=mean(mesofauna_abundance$year.c),
  sqrt_inv_n_tilda=mean(mesofauna_abundance$sqrt_inv_n_tilda)))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M4_formula,data=meso_new_data)

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(mesofauna_M4,newmods=predgrid))

#attach predictions to variables for plotting
meso_new_data <- cbind(meso_new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot data for model


#add data on sample size for each group
sample_size_label<-mesofauna_abundance%>%
  group_by(acari_collembola)%>%
  summarise(k=length(yi),study_n=n_distinct(Study_ID))%>%
  mutate(k_label=paste("k = ",k," (",study_n,")",sep = ""))

#new version of size and annual change plot
  ggplot(meso_new_data,aes(x=perc_annual_dist,y=pred,colour=acari_collembola,fill=acari_collembola))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_rep_wrap(~acari_collembola,repeat.tick.labels = TRUE)+
  geom_point(data=mesofauna_abundance,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  xlab("Change in annual precipitation (%)")+
  ylab("Change in soil fauna abundance (log response ratio)")+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)+
  geom_text(data=sample_size_label,aes(x=150,y=4,label=k_label),colour="black")

#save plot
ggsave("figures/collembola_acari.png",width = 20,height = 10,units = "cm",dpi = 300)



###############################################################################
#2 - models of heterogeneity in response of diversity of soil and litter fauna#
###############################################################################

#2.1 - taxonomic diversity

#to test models I don't want to build models with more that 4 parameters due 
#to the small dataset (k=47) and the associated risk of overparameterisation

#test all models against each other
rich_M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Study_ID/Site_ID/obsID,data=richness_complete,)
rich_M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~year.c,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~year.c+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M6<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M7<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M8<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M9<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M10<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M11<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M12<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M13<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M14<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M15<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M16<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M17<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M18<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
rich_M19<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)

#create a model selection table of these models
rich_model_sel<-model.sel(rich_M1,rich_M2,rich_M3,rich_M4,rich_M5,rich_M6,rich_M7,rich_M8,rich_M9,rich_M10,rich_M11,rich_M12,
                          rich_M13,rich_M14,rich_M15,rich_M16,rich_M18,rich_M19)

#M4 is the best model  - it just includes change in precipitation
#as well as correction for year effect and small study effect

#create dataframe for all models
rich_model_sel_df<-data.frame(rich_model_sel)
rich_model_sel_df$model<-row.names(rich_model_sel_df)

#save this model table
richness_model_sel_table<-rich_model_sel_df%>%
  select(mods,df,logLik,AICc,delta,model)%>%
  mutate(across(where(is.numeric), round, 2))%>%
  gt()
#export this to a word file
richness_model_sel_table%>%gtsave("figures/richness_selection_table.docx")

#subset to give top models
richness_model_sel_sub<-subset(rich_model_sel, delta <= 2, recalc.weights=TRUE)

#calculate model averaged coefficients
richness_coefs<-model.avg(richness_model_sel_sub)
summary_richness_coefs<-summary(richness_coefs)

#save the full averaged coefficient matrix
richness_coef_full<-data.frame(summary_richness_coefs$coefmat.full)

#tidy up this table
richness_coef_full$Variable<-row.names(richness_coef_full)
richness_coef_full_table<-richness_coef_full%>%
  relocate(Variable)%>%
  mutate(Estimate=round(Estimate,4),
         Std..Error=round(Std..Error,4),
         z.value=round(z.value,4),
         Pr...z..=round(Pr...z..,4))%>%
  rename(SE=Std..Error,
         p_value=Pr...z..)%>%
  remove_rownames()%>%
  mutate(Variable=str_replace(Variable,"perc_annual_dist","Precipitation change"),
         Variable=str_replace(Variable,"year.c","Decline effect"),
         Variable=str_replace(Variable,"intrcpt","Intercept"))%>%
  gt()

richness_coef_full_table%>%gtsave("figures/richness_coef_table.docx")

#save model formula
M4_formula<-(~perc_annual_dist)

#create new dataset for predictions
new_data_rich<-data.frame(perc_annual_dist=seq(min(richness_complete$perc_annual_dist),
                                                 max(richness_complete$perc_annual_dist),1))

#create a model matrix and remove the intercept
predgrid_rich<-model.matrix(M4_formula,data=new_data_rich)[,-1]

#predict onto the new model matrix
mypreds_rich<-data.frame(predict.rma(rich_M4,newmods=predgrid_rich))

#attach predictions to variables for plotting
new_data_rich <- cbind(new_data_rich, mypreds_rich[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#work out k and number of studies here too
richness_sample_size<-richness_complete%>%
  summarise(k=length(yi),study_n=n_distinct(Study_ID))%>%
  mutate(k_label=paste("k = ",k," (",study_n,")",sep = ""))


#plot the results of the predictions
richness_figure<-ggplot(new_data_rich,aes(perc_annual_dist,y=pred))+
  geom_line(colour="black")+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA,fill="darkgrey")+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA,fill="darkgrey")+
  geom_point(data=richness_complete,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),colour="black",alpha=0.6,shape=1)+
  theme_cowplot()+
  labs(x="Change in annual precipitation (%)",
       y="Change in soil fauna taxonomic\nrichness(log response ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)+
  geom_text(data=richness_sample_size,aes(x=30,y=1,label=k_label),colour="black")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12))
richness_figure
ggsave("figures/Figure_4.pdf",width = 12,height = 8,units = "cm",dpi = 300)



#2.2 - shannon diversity

#to test models I don't want to build models with more that 4 parameters due 
#to the small dataset (k=43) and the associated risk of overparameterisation

#test all models against each other
shannon_M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~year.c,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~year.c+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M6<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M7<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c+sqrt_inv_n_tilda,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M8<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M9<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M10<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M11<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M12<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M13<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M14<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M15<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M16<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M17<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M18<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
shannon_M19<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*exoskeleton+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)


#create a model selection table of these models
shannon_model_sel<-model.sel(shannon_M1,shannon_M2,shannon_M3,shannon_M4,shannon_M5,shannon_M6,
                             shannon_M7,shannon_M8,shannon_M9,shannon_M10,shannon_M11,shannon_M12,shannon_M13,
                             shannon_M14,shannon_M15,shannon_M16,shannon_M17,shannon_M19)

#create dataframe for all models
shannon_model_sel_df<-data.frame(shannon_model_sel)
shannon_model_sel_df$model<-row.names(shannon_model_sel)


#save this model table
shannon_model_sel_table<-shannon_model_sel_df%>%
  select(mods,df,logLik,AICc,delta,model)%>%
  mutate(across(where(is.numeric), round, 2))%>%
  gt()
#export this to a word file
shannon_model_sel_table%>%gtsave("figures/for_paper/shannon_selection_table.docx")

#subset to give top models
shannon_model_sel_sub<-subset(shannon_model_sel, delta <= 2, recalc.weights=TRUE)

#calculate model averaged coefficients
shannon_coefs<-model.avg(shannon_model_sel_sub)
summary_shannon_coefs<-summary(shannon_coefs)

#save the full averaged coefficient matrix
shannon_coef_full<-data.frame(summary_shannon_coefs$coefmat.full)

#tidy up this table
shannon_coef_full$Variable<-row.names(shannon_coef_full)
shannon_coef_full_table<-shannon_coef_full%>%
  relocate(Variable)%>%
  mutate(Estimate=round(Estimate,4),
         Std..Error=round(Std..Error,4),
         z.value=round(z.value,4),
         Pr...z..=round(Pr...z..,4))%>%
  rename(SE=Std..Error,
         p_value=Pr...z..)%>%
  remove_rownames()%>%
  mutate(Variable=str_replace(Variable,"perc_annual_dist","Precipitation change"),
         Variable=str_replace(Variable,"year.c","Decline effect"),
         Variable=str_replace(Variable,"intrcpt","Intercept"))%>%
  gt()

shannon_coef_full_table%>%gtsave("figures/for_paper/shannon_coef_table.docx")

#these results offer weak support for the hypothesis that there has been a change in impact of drought studies over time

#plot figure showing these results

#save model formula
M1_formula<-(~year.c)

#create new dataset for predictions for year
new_data_shannon_year<-data.frame(year.c=seq(min(shannon_complete$year.c),max(shannon_complete$year.c),0.1))

#create a model matrix and remove the intercept
predgrid_shannon_year<-model.matrix(M1_formula,data=new_data_shannon_year)[,-1]

#predict onto the new model matrix
mypreds_shannon_year<-data.frame(predict.rma(shannon_M1,newmods=predgrid_shannon_year))

#attach predictions to variables for plotting
new_data_shannon_year <- cbind(new_data_shannon_year, mypreds_shannon_year[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

new_data_shannon_year$study_year<-new_data_shannon_year$year.c+mean(shannon_complete$study_year)

#plot the results of the predictions
shannon_year_figure<-ggplot(new_data_shannon_year,aes(study_year,y=pred))+
  geom_line(colour="#2596be")+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA,fill="#2596be")+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA,fill="#2596be")+
  geom_point(data=shannon_complete,aes(x=study_year,y=lnrr_laj,size=1/v_lnrr_laj),colour="#2596be",alpha=0.2)+
  theme_cowplot()+
  labs(x="Year of publication",
       y="Change in soil fauna Shannon\ndiversity(log response ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)
shannon_year_figure

save_plot("figures/shannon_year_figure.png",shannon_year_figure,base_height = 10,base_width = 15,units="cm",dpi=300)


