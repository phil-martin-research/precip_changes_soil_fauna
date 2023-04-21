
# Plot the relationship between abundance and strength of disturbance for both 
# drought and precipitation increase (all) for the aridity categories 

#---------------------------------------------------
rm(list = ls())

invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

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
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)


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

#test all models against each other
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

#check to see which model is the best fit
AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11)

#copy and save the model formula
newform<-(~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1)

env_ab%>%
  group_by(Functional_group_size)%>%
  summarise(min_precip=min(perc_annual_dist),max_precip=max(perc_annual_dist))

#create dataframe with new data for predictions

new_data<-data.frame(expand.grid(perc_annual_dist = seq(-100,239,0.1),Functional_group_size=levels(as.factor(env_ab$Functional_group_size))))

#create a model matrix and remove the intercept
predgrid<-model.matrix(newform,data=new_data)

#predict onto the new model matrix
mypreds<-predict.rma(M10,newmods=predgrid)

#attach predictions to variables for plotting
new_data$pred<-mypreds$pred
new_data$ci.lb<-mypreds$ci.lb
new_data$ci.ub<-mypreds$ci.ub
new_data$pi.lb<-mypreds$pi.lb
new_data$pi.ub<-mypreds$pi.ub

#plot data for model
#first subset to remove very large effect sizes
env_ab_sub<-env_ab%>%mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
            filter(yi<2)

#remove predictions for microfauna where the change in precipitation in outside the observed values
new_data_micro<-subset(new_data,Functional_group_size=="microfauna"&perc_annual_dist<100)
new_data_no_micro<-subset(new_data,Functional_group_size!="microfauna")
#merge these together
new_data_merge<-rbind(new_data_micro,new_data_no_micro)


#new version of size and annual change plot
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
ggplot(aes(x=perc_annual_dist,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
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
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
  ggplot(aes(x=perc_annual_dist,y=(exp(pred)-1)*100,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(ci.ub)-1)*100,ymin=(exp(ci.lb)-1)*100),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(pi.ub)-1)*100,ymin=(exp(pi.lb)-1)*100),colour=NA)+
  facet_wrap(~Functional_group_size,scales = "free")+
  geom_point(data=env_ab_sub,aes(x=perc_annual_dist,y=(exp(yi)-1)*100,size=1/vi),alpha=0.25)+
  xlab("change in annual precipitation (%)")+
  ylab("change in soil fauna abundance (%)")+
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figures/abundance_precip_size.png",width = 18,height = 8,units = "cm",dpi = 300)

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


#create values over which we want to generate predictions
perc_annual_dist<-values(europe_perc_change_forest)#percentage annual change in precipitation in European forests



Functional_group_sizemacrofauna<-seq(0,0,length=length(perc_annual_dist))
Functional_group_sizemesofauna<-seq(1,1,length=length(perc_annual_dist))
Functional_group_sizemicrofauna<-seq(0,0,length=length(perc_annual_dist))
perc_annual_dist2<-values(europe_perc_change_forest)^2
perc_meso<-values(europe_perc_change_forest)*1
perc_micro<-0
perc2_meso<-1
perc2_micro<-0
newmods<-cbind(perc_annual_dist,
               Functional_group_sizemacrofauna,Functional_group_sizemesofauna,Functional_group_sizemicrofauna,
               perc_annual_dist2,perc_meso,perc_micro,perc2_meso,perc2_micro)




model.matrix(M10)

predict(M10,newmods=newmods)




