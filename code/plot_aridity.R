
# Plot the relationship between abundance and strength of disturbance for both 
# drought and precipitation increase (all) for the aridity categories 

#---------------------------------------------------
rm(list = ls())

library(ggpubr)
library(cowplot)
library(readxl)
library(ggplot2)
library(metafor)
library(dplyr)
library(tidyverse)
library(readr)
library(forcats)
library(raster)
library(maptools)

install.packages("vctrs")

Env_df <- read_csv("data/ALL_env_df.csv")

#---------------------------------------------------
#adapt dataset for analysis
#---------------------------------------------------

#remove anomalies - outlying disturbance strengths 
Env_df$perc_annual_dist <- ifelse(Env_df$perc_annual_dist<500, Env_df$perc_annual_dist,  NA)
Env_df$perc_during_dist <- ifelse(Env_df$perc_during_dist<500, Env_df$perc_during_dist,  NA)

#change variances from 0 to NA
Env_df$control_var <- ifelse(Env_df$control_var !=0, Env_df$control_var,  NA)
Env_df$dist_var <- ifelse(Env_df$dist_var !=0, Env_df$dist_var, NA)

# change percentage to be continuous from decrease to increase 
Env_df$perc_annual_dist <- ifelse(Env_df$perc_annual_dist>100, Env_df$perc_annual_dist-100, Env_df$perc_annual_dist*(-1))
Env_df$perc_during_dist <- ifelse(Env_df$perc_during_dist>100, Env_df$perc_during_dist-100, Env_df$perc_during_dist*(-1))

#convert SE to SD
Env_df$control_SD<-ifelse(Env_df$var_type=="SE", Env_df$control_var*sqrt(Env_df$control_n), Env_df$control_var)
Env_df$dist_SD<-ifelse(Env_df$var_type=="SE",Env_df$dist_var*sqrt(Env_df$dist_n), Env_df$dist_var)


#---------------------------------------------------
#1. data exploration and effect sizes
#---------------------------------------------------
#log response ratio
Env_df_rr <- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                    sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = Env_df)

#impute missing variances with median variance 
Env_df_rr$vi<-ifelse(is.na(Env_df_rr$vi),median(Env_df_rr$vi, na.rm = TRUE), Env_df_rr $vi)

#subset Env_df_rr dataset to set variables of interest
env_ab_rr <- Env_df_rr %>%filter(broad_outcome == 'abundance')


#complete cases of variables and yi and filter cooks 
complete_env_ab_rr <- env_ab_rr[complete.cases(env_ab_rr[ , c(33:35,79)]),]

complete_env_ab_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=complete_env_ab_rr)
cooks<- cooks.distance(complete_env_ab_m0)
env_ab <- complete_env_ab_rr %>%cbind(cooks) %>%filter(cooks < 3.0*mean(cooks))

# models 

env_ab_m2<-rma.mv(yi,vi,mods = ~aridity,random=~1|Site_ID/Study_ID,data=env_ab)

env_df_m2_df <- as.data.frame(predict(env_ab_m2, newmods=c(0:3), addx=TRUE))
env_plot_m2 <- ggplot()+
  geom_point(data = env_ab, aes(x = aridity , y = yi, size = 1/vi), alpha =0.5)+
  geom_ribbon(data = env_df_m2_df, 
              aes(x = X.aridity, y = pred, ymin = ci.lb, ymax = ci.ub), alpha = 0.2)+
  geom_line(data = env_df_m2_df, 
            aes(x = X.aridity, y = pred), col= 'red')+
  ggtitle("aridity")+
  theme_cowplot()


env_ab$cat <- cut(env_ab$aridity,
                  breaks=c(0,.05, .2, .5, .65,3),
                  labels=c('hyper_arid', 'arid', 'semi_arid', 'dry_sub_humid','humid' ))


#set variables of interest 
#no hyper_arid or arid regions -> 
#hyper_arid<- env_ab %>% filter(broad_outcome == 'abundance')%>%filter(cat == 'hyper_arid')
#arid <- env_ab %>%filter(broad_outcome == 'abundance')%>%filter(cat == 'arid')
semi_arid <- env_ab %>% filter(broad_outcome == 'abundance')%>%filter(cat == 'semi_arid')
dry_sub_humid <- env_ab %>% filter(broad_outcome == 'abundance')%>%filter(cat == 'dry_sub_humid')
humid <- env_ab %>% filter(broad_outcome == 'abundance')%>%filter(cat == 'humid')

####################################################
#1 - data exploration and effect sizes##############
####################################################
#log response ratio
semi_arid_rr <- escalc(m2i = control_av, m1i = disturbance_av,n2i = control_n, n1i = dist_n,
                       sd2i = control_SD, sd1i = dist_SD,  measure = "ROM",  data = semi_arid)
dry_sub_humid_rr <- escalc(m2i = control_av, m1i = disturbance_av,n2i = control_n, n1i = dist_n,
                           sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = dry_sub_humid )
humid_rr <- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                   sd2i = control_SD, sd1i = dist_SD, measure = "ROM",  data = humid )

#impute variation 
semi_arid_rr$vi<-ifelse(is.na(semi_arid_rr$vi),median(semi_arid_rr$vi, na.rm = TRUE), semi_arid_rr $vi)
dry_sub_humid_rr$vi<-ifelse(is.na(dry_sub_humid_rr$vi), median(dry_sub_humid_rr$vi, na.rm = TRUE), dry_sub_humid_rr$vi)
humid_rr$vi<-ifelse(is.na(humid_rr$vi), median(humid_rr$vi, na.rm = TRUE),humid_rr$vi)

#complete cases for yi
semi_arid_rr <- semi_arid_rr[complete.cases(semi_arid_rr[ ,79 ]),]
dry_sub_humid_rr <- dry_sub_humid_rr[complete.cases(dry_sub_humid_rr[ ,79 ]),]
humid_rr <- humid_rr[complete.cases(humid_rr[ ,79 ]),]

# models 
semi_arid_m0 <-rma.mv(yi,vi, random=~1|Study_ID/Site_ID,data=semi_arid_rr)
semi_arid_m1 <-rma.mv(yi,vi,mods = ~perc_annual_dist, random=~1|Study_ID/Site_ID,data=semi_arid_rr)
semi_arid_m2 <-rma.mv(yi,vi,mods = ~time_after_dist_start, random=~1|Study_ID/Site_ID,data=semi_arid_rr)
dry_sub_humid_m0 <-rma.mv(yi,vi, random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
dry_sub_humid_m1 <-rma.mv(yi,vi,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
dry_sub_humid_m2 <-rma.mv(yi,vi,mods = ~time_after_dist_start,random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
humid_m0 <-rma.mv(yi,vi, random=~1|Study_ID/Site_ID,data=humid_rr)
humid_m1 <-rma.mv(yi,vi,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID,data=humid_rr)
humid_m2 <-rma.mv(yi,vi,mods = ~time_after_dist_start,random=~1|Study_ID/Site_ID,data=humid_rr)



dry_sub_humid_m2<-rma.mv(yi,vi,mods = ~disturbance_type-1, random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
semi_arid_m2 <-rma.mv(yi,vi,mods = ~disturbance_type-1,random=~1|Study_ID/Site_ID,data=semi_arid_rr)
humid_m2 <-rma.mv(yi,vi,mods = ~disturbance_type-1,random=~1|Study_ID/Site_ID,data=humid_rr)

dry_sub_humid_m3<-rma.mv(yi,vi,mods = ~perc_annual_dist, random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
semi_arid_m3 <-rma.mv(yi,vi,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID,data=semi_arid_rr)
humid_m3 <-rma.mv(yi,vi,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID,data=humid_rr)

dry_sub_humid_m4<-rma.mv(yi,vi,mods = ~disturbance_type, random=~1|Study_ID/Site_ID,data=dry_sub_humid_rr)
semi_arid_m4 <-rma.mv(yi,vi,mods = ~disturbance_type,random=~1|Study_ID/Site_ID,data=semi_arid_rr)
humid_m4 <-rma.mv(yi,vi,mods = ~disturbance_type,random=~1|Study_ID/Site_ID,data=humid_rr)


pred_dry_sub_humid_m1_df <- as.data.frame(predict(dry_sub_humid_m1, newmods=c(-100:50), addx=TRUE))
dry_sub_humid1 <- ggplot()+
  geom_point(data = dry_sub_humid_rr, aes(x = perc_annual_dist , y = yi, size = 1/vi), alpha =0.5)+
  geom_ribbon(data = pred_dry_sub_humid_m1_df, 
              aes(x = X.perc_annual_dist, y = pred, ymin = ci.lb, ymax = ci.ub), alpha = 0.2)+
  geom_line(data = pred_dry_sub_humid_m1_df, 
            aes(x = X.perc_annual_dist, y = pred), col= 'red')+
  ggtitle("dry_sub_humid")+
  ylim(-5,3)+
  theme_cowplot()

pred_semi_arid_m1_df <- as.data.frame(predict(semi_arid_m1, newmods=c(-30:3), addx=TRUE))
semi_arid1 <- ggplot()+
  geom_point(data = semi_arid_rr, aes(x = perc_annual_dist , y = yi, size = 1/vi), alpha =0.5)+
  geom_ribbon(data = pred_semi_arid_m1_df, 
              aes(x = X.perc_annual_dist, y = pred, ymin = ci.lb, ymax = ci.ub), alpha = 0.2)+
  geom_line(data = pred_semi_arid_m1_df, 
            aes(x = X.perc_annual_dist, y = pred), col= 'red')+
  ylim(-5,3)+
  ggtitle("semi_arid")+
  theme_cowplot()

pred_humid_m1_df <- as.data.frame(predict(humid_m1, newmods=c(-100:250), addx=TRUE))
humid1 <- ggplot()+
  geom_point(data = humid_rr, aes(x = perc_annual_dist , y = yi, size = 1/vi), alpha =0.5)+
  geom_ribbon(data = pred_humid_m1_df, 
              aes(x = X.perc_annual_dist, y = pred, ymin = ci.lb, ymax = ci.ub), alpha = 0.2)+
  geom_line(data = pred_humid_m1_df, 
            aes(x = X.perc_annual_dist, y = pred), col= 'red')+
  ylim(-5,3)+
  ggtitle("humid")+
  theme_cowplot()

ggarrange(semi_arid1 , dry_sub_humid1  + rremove("ylab"),humid1+ rremove("ylab"),
          nrow = 1)







aridity_index_df <- aridity_continuous_df %>% 
  mutate(category = case_when(
    is.infinite(layer) ~ 'Humid',
    layer >= 0.65 ~ 'Humid',
    layer >= 0.5 & layer < 0.65 ~ 'Dry sub-humid',
    layer >= 0.2 & layer < 0.5 ~ 'Semi-arid',
    layer >= 0.05 & layer < 0.2 ~ 'Arid',
    layer < 0.05 ~ 'Hyper-arid'
  )) %>% 
  # Convert to ordered factor
  mutate(category = factor(category,
                           levels = c('Hyper-arid', 'Arid', 'Semi-arid',
                                      'Dry sub-humid', 'Humid'),
                           ordered = TRUE))



########################################################################
#Phil's models##########################################################
########################################################################
names(env_ab)

ggplot(env_ab,aes(perc_annual_dist,yi))+
  geom_point()

M0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=env_ab)
M1<-rma.mv(yi,vi,mods = ~aridity,random=~1|Site_ID/Study_ID,data=env_ab)
M2<-rma.mv(yi,vi,mods = ~perc_annual_dist ,random=~1|Site_ID/Study_ID,data=env_ab)
M3<-rma.mv(yi,vi,mods = ~Functional_group_size ,random=~1|Site_ID/Study_ID,data=env_ab)
M4<-rma.mv(yi,vi,mods = ~aridity+perc_annual_dist,random=~1|Site_ID/Study_ID,data=env_ab)
M5<-rma.mv(yi,vi,mods = ~aridity*perc_annual_dist,random=~1|Site_ID/Study_ID,data=env_ab)
M6<-rma.mv(yi,vi,mods = ~aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=env_ab)
M7<-rma.mv(yi,vi,mods = ~aridity*Functional_group_size,random=~1|Site_ID/Study_ID,data=env_ab)
M8<-rma.mv(yi,vi,mods = ~perc_annual_dist+Functional_group_size,random=~1|Site_ID/Study_ID,data=env_ab)
M9<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=env_ab)
M10<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=env_ab)
M11<-rma.mv(yi,vi,mods = ~perc_annual_dist*aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=env_ab)

AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11)

summary(env_ab$aridity)
summary(env_ab$perc_annual_dist)
arid_prec_df<-expand.grid(aridity=seq(0.32,2.16,0.01),perc_annual_dist=seq(-100,239,1))

M9_preds<-data.frame(predict(M10,addx=TRUE))

head(M9_preds)

unique(env_ab$Functional_group_size)

M9_preds$Functional_group_size<-ifelse(M9_preds$X.Functional_group_sizemacrofauna==1,"macrofauna",NA)
M9_preds$Functional_group_size<-ifelse(M9_preds$X.Functional_group_sizemesofauna==1,"mesofauna",M9_preds$Functional_group_size)
M9_preds$Functional_group_size<-ifelse(M9_preds$X.Functional_group_sizemicrofauna==1,"microfauna",M9_preds$Functional_group_size)

env_ab_sub<-env_ab%>%mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
            filter(yi<2)


#new version of size and annual change plot
M9_preds%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
ggplot(aes(x=X.perc_annual_dist,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_wrap(~Functional_group_size,scales = "free")+
  geom_point(data=env_ab_sub,aes(x=perc_annual_dist,y=yi,size=1/vi),alpha=0.25)+
  xlab("change in annual precipitation (%)")+
  ylab("change in soil fauna abundance (log ratio)")+
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

#same plot but with percentage change on the abundance axis
M9_preds%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
  ggplot(aes(x=X.perc_annual_dist,y=(exp(pred)-1)*100,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(ci.ub)-1)*100,ymin=(exp(ci.lb)-1)*100),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(pi.ub)-1)*100,ymin=(exp(pi.lb)-1)*100),colour=NA)+
  facet_wrap(~Functional_group_size,scales = "free")+
  geom_point(data=env_ab_sub,aes(x=perc_annual_dist,y=(exp(yi)-1)*100,size=1/vi),alpha=0.25)+
  xlab("change in annual precipitation (%)")+
  ylab("change in soil fauna abundance (log ratio)")+
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#look at how much data there is for different groups
env_ab%>%group_by(Taxa_subclass)%>%
  summarise(total=length(yi))%>%
  print(n=26)


#subset just for collembola
env_ab_coll<-env_ab%>%
  filter(Taxa_subclass=="Collembola")


ggplot(env_ab_coll,aes(perc_annual_dist,yi))+
  geom_point()+
  geom_smooth()

#subset just for acari
env_ab_acari<-env_ab%>%
  filter(Taxa_subclass=="Acari")


ggplot(env_ab_acari,aes(perc_annual_dist,yi))+
  geom_point()


#read in spatial data

precip_2000<-raster("data/spatial_data/stacked-mmyr-abs-annual_pr_rcp45_ens_1985-2015_nexgddp.tif")

precip_2070<-raster("data/spatial_data/stacked-mmyr-abs-annual_pr_rcp45_ens_2065-2095_nexgddp.tif")


plot((precip_2070-precip_2000)/precip_2000)

plot(wrld_simpl)
