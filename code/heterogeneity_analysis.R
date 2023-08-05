# script to explore heterogeneity in the effects of precipitation changes on soil and litter fauna 

rm(list = ls())

pacman::p_load(tidyverse,metafor,cowplot,orchaRd,ggbeeswarm,tidyr,ggthemes,sp,broom,lemon,MuMIn,glmulti,PerformanceAnalytics,GGally,gt)

#read in .csv files with soil fauna data
abundance<- read_csv("data/abundance_data.csv")
richness<- read_csv("data/richness_data.csv")
shannon<- read_csv("data/shannon_data.csv")


###############################################################################
#1 - data tidying##############################################################
###############################################################################

abundance%>%
  group_by(Validity,perc_annual_dist,Functional_group_size.y,aridity,randomized)%>%
  summarise(size_count=length(above_below))%>%
  print(n=100)

#complete cases for the variable about percentage change in precipitation and body size
abundance_complete<-abundance[complete.cases(abundance$perc_annual_dist,abundance$Functional_group_size.y,abundance$above_below,
                                             abundance$exoskeleton),]
#we lose 21 comparisons for abundance

richness_complete<-richness[complete.cases(richness$perc_annual_dist,
                                           richness$Functional_group_size.y,
                                           richness$above_below,
                                           richness$exoskeleton),]
#we lose 5 comparisons for richness

shannon_complete<-shannon[complete.cases(shannon$perc_annual_dist,shannon$Functional_group_size.y,shannon$above_below,
                                         shannon$exoskeleton),]
#we lose 2 comparisons for shannon wiener


###############################################################################
#1 - models of heterogeneity in response of abundance of soil and litter fauna#
###############################################################################

# create a unit-level random effect to model residual variance in metafor
abundance_complete$obsID <- 1:nrow(abundance_complete)
# mean-centering year of publication to help with interpretation
abundance_complete$year.c <- as.vector(scale(abundance_complete$study_year, scale = F))

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

eval(metafor:::.MuMIn)

#create a model selection table of these models
abun_model_sel<-model.sel(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,M16)

#M12 is the best model that includes an interaction between functional group size and change in precipitation
#as well as correction for year effect and small study effect

#create dataframe for all models
abundance_model_sel_df<-data.frame(abun_model_sel)
abundance_model_sel_df$model<-row.names(abundance_model_sel_df)

#calculate R2 for all models using a loop
abundance_model_list<-list(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,M16)
abundance_model_list_names<-c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M14","M15","M16")
abundance_model_R2<-NULL
for (i in 1:length(abundance_model_list)){
  abundance_model_R2_temp<-data.frame(model=abundance_model_list_names[i],
                                      R2=(sum(M0$sigma2) - sum(abundance_model_list[[i]]$sigma2)) / sum(M0$sigma2))
  abundance_model_R2<-rbind(abundance_model_R2,abundance_model_R2_temp)
}

#merge R2 with other model descriptors
abundance_model_sel_df2<-abundance_model_sel_df%>%
  left_join(abundance_model_R2,"model")%>%
  arrange(desc(weight))

#save this model table
abundance_model_sel_table<-abundance_model_sel_df2%>%
  mutate(across(where(is.numeric), round, 3))%>%
  gt()
#expord this to a word file
abundance_model_sel_table%>%gtsave("figures/for_paper/abundance_selection_table.docx")


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
  ylab("Change in soil fauna abundance (log ratio)")+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_fill_manual(values = c("#02475f","#c3386b","#e0b500"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)
ggsave("figures/for_paper/abundance_precip_size.png",width = 20,height = 10,units = "cm",dpi = 300)


###############################################################################
#2 - models of heterogeneity in response of diversity of soil and litter fauna#
###############################################################################

#2.1 - taxonomic diversity

# create a unit-level random effect to model residual variance in metafor
richness_complete$obsID <- 1:nrow(richness_complete)
# mean-centering year of publication to help with interpretation
richness_complete$year.c <- as.vector(scale(richness_complete$study_year, scale = F))

#first a saturated model including all potential predictors
M_sat_rich<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+
                                              Functional_group_size.y+
                                              above_below+
                                              exoskeleton+
                                              year.c+
                                              sqrt_inv_n_tilda-1,
                                              random=~1|Study_ID/Site_ID/obsID,data=richness_complete)


#to test models I don't want to build models with more that 4 parameters due 
#to the small dataset (k=43) and the associated risk of overparameterisation


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
                          rich_M13,rich_M14,rich_M15,rich_M16,rich_M17,rich_M18,rich_M19)

#M4 is the best model  - it just includes change in precipitation
#as well as correction for year effect and small study effect

#create dataframe for all models
rich_model_sel_df<-data.frame(rich_model_sel)
rich_model_sel_df$model<-row.names(rich_model_sel_df)

#calculate R2 for all models using a loop
rich_model_list<-list(rich_M1,rich_M2,rich_M3,rich_M4,rich_M5,rich_M6,rich_M7,rich_M8,rich_M9,rich_M10,rich_M11,rich_M12,
                           rich_M13,rich_M14,rich_M15,rich_M16,rich_M17,rich_M18,rich_M19)
rich_model_list_names<-c("rich_M1","rich_M2","rich_M3","rich_M4","rich_M5","rich_M6","rich_M7","rich_M8","rich_M9","rich_M10","rich_M11","rich_M12",
                              "rich_M13","rich_M14","rich_M15","rich_M16","rich_M17","rich_M18","rich_M19")
rich_model_R2<-NULL
for (i in 1:length(rich_model_list)){
  rich_model_R2_temp<-data.frame(model=rich_model_list_names[i],
                                      R2=(sum(rich_M0$sigma2) - sum(rich_model_list[[i]]$sigma2)) / sum(rich_M0$sigma2))
  rich_model_R2<-rbind(rich_model_R2,rich_model_R2_temp)
}

#set negative R2 values to zero
rich_model_R2$R2<-ifelse(rich_model_R2$R2<0,0,rich_model_R2$R2)

#merge R2 with other model descriptors
rich_model_sel_df2<-rich_model_sel_df%>%
  left_join(rich_model_R2,"model")%>%
  arrange(desc(weight))

#save this model table
rich_model_sel_table<-rich_model_sel_df2%>%
  mutate(across(where(is.numeric), round, 3))%>%
  gt()
#export this to a word file
rich_model_sel_table%>%gtsave("figures/for_paper/richness_selection_table.docx")


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


#plot the results of the predictions
richness_figure<-ggplot(new_data_rich,aes(perc_annual_dist,y=pred))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  geom_point(data=richness_complete,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  theme_cowplot()+
  labs(x="Change in annual precipitation (%)",
       y="Change in soil fauna taxonomic richness\n(log response ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)+
  ylim(-2,1)
richness_figure
ggsave("figures/for_paper/richness_precip.png",width = 12,height = 8,units = "cm",dpi = 300)



#2.1 - shannon diversity

# create a unit-level random effect to model residual variance in metafor
shannon_complete$obsID <- 1:nrow(shannon_complete)
# mean-centering year of publication to help with interpretation
shannon_complete$year.c <- as.vector(scale(shannon_complete$study_year, scale = F))

#first a saturated model including all potential predictors
M_sat_shannon<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+
                     Functional_group_size.y+
                     above_below+
                     exoskeleton+
                     year.c+
                     sqrt_inv_n_tilda-1,
                   random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)


#to test models I don't want to build models with more that 4 parameters due 
#to the small dataset (k=38) and the associated risk of overparameterisation


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


#check to see which model is the most parsimonious 
AIC.rma(shannon_M0,shannon_M1,shannon_M2,shannon_M3,shannon_M4,shannon_M5,shannon_M6,
        shannon_M7,shannon_M8,shannon_M9,shannon_M10,shannon_M11,shannon_M12,shannon_M13,
        shannon_M14,shannon_M15,shannon_M16,shannon_M17,shannon_M19)

#these results support the hypothesis that there has been a change in impact of drought studies over time
#however ecologically, it's not very interesting!

#######################################################
#mapping analysis######################################
#######################################################

#read in precipitation data
precip_2000<-raster("data/spatial_data/stacked-mmyr-abs-annual_pr_rcp45_ens_1985-2015_nexgddp.tif")
precip_2070<-raster("data/spatial_data/stacked-mmyr-abs-annual_pr_rcp45_ens_2065-2095_nexgddp.tif")
#import forest data
dec_broad<-raster("data/spatial_data/land_cover/deciduous_broadleaf.tif")
ever_broad<-raster("data/spatial_data/land_cover/evergreen_broadleaf.tif")
ever_needle<-raster("data/spatial_data/land_cover/evergreen_deciduous_needleleaf.tif")
mixed_trees<-raster("data/spatial_data/land_cover/mixed_other_trees.tif")
#read in world shapefile
data(wrld_simpl)
#create SpatialPolygon the size of Europe
eur_crop<-as(extent(-10,40,35,70),"SpatialPolygons")
#set coordinate system to be the same for Europe polygon and precipitation data
crs(eur_crop)<-crs(precip_2000)
#crop precipitation data to just Europe
precip_2000_crop<-crop(precip_2000,eur_crop)
precip_2070_crop<-crop(precip_2070,eur_crop)
#calculate percentage change in precipitation in Europe between 2000 and 2070
perc_change_crop<-(((precip_2070_crop-precip_2000_crop)/precip_2000_crop)*100)
#mask percentage change in precipitation so that there is only data for terrestrial areas
mask_perc_change_crop<-mask(perc_change_crop,wrld_simpl)
#crop forest data to just Europe
dec_broad_crop<-crop(dec_broad,eur_crop)
ever_broad_crop<-crop(ever_broad,eur_crop)
ever_needle_crop<-crop(ever_needle,eur_crop)
mixed_trees_crop<-crop(mixed_trees,eur_crop)
#add together all data on trees to give total cover
all_trees_crop<-dec_broad_crop+ever_broad_crop+ever_needle_crop+mixed_trees_crop
#reclassify so that forest is defined as areas where tree cover is >40%
m<-c(-Inf,40,NA,40,102,1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
forest_40<-reclassify(all_trees_crop, rclmat)
#resample precipitation to have the same resolution as forest cover data
precip_resample<-terra::resample(mask_perc_change_crop,forest_40)
#crop precipitation change data so that it only represents changes in forest areas
forest_precip<-precip_resample*forest_40

#turn this into a dataframe
forest_precip_df<-as.data.frame(forest_precip,xy=TRUE)

#subset the dataframe to just Europe to test model
forest_precip_df_Spain<-subset(forest_precip_df,x>-10&x<30&y<80&y>36)
forest_precip_df_Spain<-subset(forest_precip_df_Spain,!is.na(layer))

ggplot(forest_precip_df_Spain,aes(x,y,fill=layer))+
geom_tile()

#make predictions for the map
#copy and save the model formula
newform<-(~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1)

#take values from map over which we want to generate predictions
perc_annual_dist_raster<-forest_precip_df_Spain$layer

#create dataframe with new data for predictions
new_data_map<-data.frame(expand.grid(perc_annual_dist = perc_annual_dist_raster,Functional_group_size=levels(as.factor(abundance$Functional_group_size))))

#create a model matrix and remove the intercept
predgrid_map<-model.matrix(newform,data=new_data_map)

#predict onto the new model matrix
mypreds_map<-predict.rma(M10,newmods=predgrid_map)

#attach predictions to variables for plotting
new_data_map$pred<-mypreds_map$pred
new_data_map$ci.lb<-mypreds_map$ci.lb
new_data_map$ci.ub<-mypreds_map$ci.ub
new_data_map$pi.lb<-mypreds_map$pi.lb
new_data_map$pi.ub<-mypreds_map$pi.ub

#attach latitude and longitude
new_data_map$lat<-forest_precip_df_Spain$y
new_data_map$long<-forest_precip_df_Spain$x

ggplot(new_data_map,aes(long,lat,fill=(exp(pred)-1)*100))+
  geom_tile()+
  scale_fill_gradient2()+
  facet_wrap(~Functional_group_size)

ggplot(forest_precip_df,aes(x,y,fill=layer))+
  geom_tile()


