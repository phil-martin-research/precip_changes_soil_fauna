#code for quick analysis of precipitation change data

#load packages
library(orchaRd)
library(tidyverse)
library(cowplot)
library(metafor)

#load in data
precip_data<-read.csv("data/study_data_for_analysis_07_06.csv")

#remove anomalies 
precip_data$control_var <- ifelse(precip_data$control_var!=0,   precip_data$control_var, NA)
precip_data$dist_var <- ifelse(precip_data$dist_var!=0, precip_data$dist_var, NA)

#convert standard errors and confidence intervals to SD

precip_data$control_var_new<-ifelse(precip_data$var_type=="SE",
                                    precip_data$control_var*sqrt(precip_data$control_n),
                                    precip_data$control_var)

precip_data$dist_var_new<-ifelse(precip_data$var_type=="SE",
                                    precip_data$dist_var*sqrt(precip_data$dist_n),
                                    precip_data$dist_var)

#calculate coefficient of variation for each outcome and disturbance grouping
cv_values<-precip_data%>%
  group_by(broad_outcome,disturbance_type)%>%
  summarise(cv_dist=mean(dist_var_new/disturbance_av,na.rm=TRUE),
            cv_control=mean(control_var_new/control_av,na.rm=TRUE))

#replace missing values with median based on coefficient of variation
# abundance and drought - disturbed
precip_data$dist_var_new_imputed<-ifelse(is.na(precip_data$dist_var_new)&
                                           precip_data$broad_outcome=="abundance"&
                                           precip_data$disturbance_type=="drought",
                                           0.967*precip_data$disturbance_av,
                                           precip_data$dist_var_new)

#abundance and precipitation increases - disturbed
precip_data$dist_var_new_imputed<-ifelse(is.na(precip_data$dist_var_new_imputed)&
                                           precip_data$broad_outcome=="abundance"&
                                           precip_data$disturbance_type=="precip_inc",
                                           0.858*precip_data$disturbance_av,
                                           precip_data$dist_var_new_imputed)
# abundance and drought - control
precip_data$control_var_new_imputed<-ifelse(is.na(precip_data$control_var_new)&
                                           precip_data$broad_outcome=="abundance"&
                                           precip_data$disturbance_type=="drought",
                                           0.723*precip_data$control_av,
                                           precip_data$dist_var_new)

#abundance and precipitation increases - disturbed
precip_data$control_var_new_imputed<-ifelse(is.na(precip_data$control_var_new)&
                                            precip_data$broad_outcome=="abundance"&
                                            precip_data$disturbance_type=="precip_inc",
                                            1.38*precip_data$control_av,
                                            precip_data$dist_var_new)

#replace mean values of zero with values of 0.1
precip_data$control_av_new<-ifelse(precip_data$control_av==0,0.1,precip_data$control_av)
precip_data$disturbance_av_new<-ifelse(precip_data$disturbance_av==0,0.1,precip_data$disturbance_av)


#######################################
#meta-analysis#########################
#######################################

#subset abundance data
abun<-subset(precip_data,broad_outcome=="abundance")
abun_drought<-subset(precip_data,broad_outcome=="abundance"&disturbance_type=="drought")
abun_precip<-subset(precip_data,broad_outcome=="abundance"&disturbance_type=="precip_inc")

#calculate the effect size
es_ab<-escalc("ROM",m2i=control_av_new,m1i=disturbance_av_new,
                      sd2i=control_var_new_imputed,sd1i=dist_var_new_imputed,
                      n2i=control_n,n1i=dist_n,data=abun)

es_ab_filter<-es_ab%>%filter(!is.na(yi)|!is.na(vi))

es_ab_drought<-escalc("ROM",m2i=control_av,m1i=disturbance_av,
       sd2i=control_var_new_imputed,sd1i=dist_var_new_imputed,
       n2i=control_n,n1i=dist_n,data=abun_drought)

es_ab_precip<-escalc("ROM",m2i=control_av,m1i=disturbance_av,
                          sd2i=control_var_new_imputed,sd1i=dist_var_new_imputed,
                          n2i=control_n,n1i=dist_n,data=abun_precip)

#run model
abun_model<-rma.mv(yi,vi,mods = ~disturbance_type-1,random=~1|Study_ID/Site_ID,data=es_ab)


metafor::cooks.distance.rma.mv(abun_model)

cooks_ab_model<-cooks.distance(abun_model)
cooks_ab <- es_ab %>%cbind(cooks_ab_model)%>%filter(cooks_ab_model < 3.0*mean(cooks_ab_model))

?cooks.distance



ab_d_model_cooks<-rma.mv(yi,vi,random=~1|Study_ID/Site_ID,data=cooks_d_ab)



hist(es_ab_drought$yi)

#plot model results
orchard_plot(abun_model,mod="disturbance_type",group="Study_ID",data=es_precip_data,xlab="log response ratio")
