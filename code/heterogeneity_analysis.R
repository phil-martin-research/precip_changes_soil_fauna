# script to explore heterogeneity in the effects of precipitation changes on soil and litter fauan 

rm(list = ls())

library(tidyverse)
library(cowplot)
library(metafor)
library(orchaRd)
library(ggbeeswarm)
library(tidyr)


#read in .csv files with soil fauna data
abundance<- read_csv("data/abundance_data.csv")
diversity<- read_csv("data/diversity_data.csv")

###############################################################################
#1 - imputation of missing values for precipitation change######################
###############################################################################

sum(is.na(abundance$perc_annual_dist))
#six values are missing for the magnitude of precipitation change, so we can impute the median value
#for studies of reduction and increases in precipitation

abundance%>%
  group_by(disturbance_type)%>%
  summarise(med_mag=median(perc_annual_dist,na.rm=TRUE))
#impute these values
abundance<-abundance%>%
  mutate(perc_annual_dist_type=if_else(!is.na(perc_annual_dist),"not_imputed","NA"),
         perc_annual_dist=if_else(is.na(perc_annual_dist)&disturbance_type=="drought",-34.2,perc_annual_dist),
         perc_annual_dist_type=if_else(is.na(perc_annual_dist)&disturbance_type=="drought","imputed",perc_annual_dist_type),
         perc_annual_dist=if_else(is.na(perc_annual_dist)&disturbance_type=="precip_inc",100,perc_annual_dist),
         perc_annual_dist_type=if_else(is.na(perc_annual_dist)&disturbance_type=="precip_inc","imputed",perc_annual_dist_type),
         )


#complete cases for the variable about percentage change in precipitation
abundance_complete<-abundance[complete.cases(abundance$perc_annual_dist),]
#we lose 6 comparisons for abundance
diversity_complete<-diversity[complete.cases(diversity$perc_annual_dist),]
#we lose 4 comparisons for diversity

###############################################################################
#1 - models of heterogeneity in response of abundance of soil and litter fauna#
###############################################################################

#first a saturated model including all potential predictors
M_sat_abun<-rma.mv(yi,vi,mods = ~aridity+perc_annual_dist+I(perc_annual_dist^2)+Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_complete)
#use cooks distance to identify influential points
cooks_sat_abun<-cooks.distance(M_sat_abun)
#filter out high cooks distances for saturated model and then run all models
abundance_filtered<- abundance_complete %>%cbind(cooks_sat_abun) %>%filter(cooks_sat_abun < 3.0*mean(cooks_sat_abun,na.rm=TRUE))
#this removes 5 comparisons

#test all models against each other
M0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M1<-rma.mv(yi,vi,mods = ~aridity,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M2<-rma.mv(yi,vi,mods = ~perc_annual_dist ,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M3<-rma.mv(yi,vi,mods = ~Functional_group_size ,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M4<-rma.mv(yi,vi,mods = ~aridity+perc_annual_dist,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M5<-rma.mv(yi,vi,mods = ~aridity*perc_annual_dist,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M6<-rma.mv(yi,vi,mods = ~aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M7<-rma.mv(yi,vi,mods = ~aridity*Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M8<-rma.mv(yi,vi,mods = ~perc_annual_dist+Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M9<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M10<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=abundance_filtered)
M11<-rma.mv(yi,vi,mods = ~perc_annual_dist*aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_filtered)

#check to see which model is the best fit
AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11)

#copy and save the model formula
M10_formula<-(~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1)

#create dataframe with new data for predictions
new_data<-data.frame(expand.grid(perc_annual_dist = seq(-100,239,0.1),Functional_group_size=levels(as.factor(abundance_filtered$Functional_group_size))))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M10_formula,data=new_data)

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(M10,newmods=predgrid))

#attach predictions to variables for plotting
new_data <- cbind(new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot data for model
#first subset to remove very large effect sizes
abundance_sub<-abundance_filtered%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
            filter(yi<2)

#remove predictions for microfauna and macrofauna where the changes in precipitation are outside the observed values
new_data_micro<-subset(new_data,Functional_group_size=="microfauna"&perc_annual_dist<100)
new_data_meso<-subset(new_data,Functional_group_size=="mesofauna")
new_data_macro<-subset(new_data,Functional_group_size=="macrofauna"&perc_annual_dist<100)

#merge these together
new_data_merge<-rbind(new_data_micro,new_data_meso,new_data_macro)


#new version of size and annual change plot
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
ggplot(aes(x=perc_annual_dist,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_wrap(~Functional_group_size,scales = "free")+
  geom_point(data=abundance_sub,aes(x=perc_annual_dist,y=yi,size=1/vi),alpha=0.25)+
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
  geom_point(data=abundance_sub,aes(x=perc_annual_dist,y=(exp(yi)-1)*100,size=1/vi),alpha=0.25)+
  xlab("change in annual precipitation (%)")+
  ylab("change in soil fauna abundance (%)")+
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("figures/for_paper/abundance_precip_size.png",width = 20,height = 10,units = "cm",dpi = 300)


#check for impact of publication bias on this
# create a unit-level random effect to model residual variance in metafor
abundance_filtered$obsID <- 1:nrow(abundance_filtered)
# mean-centering year of publication to help with interpretation
abundance_filtered$year.c <- as.vector(scale(abundance_filtered$study_year, scale = F))

#run an all-in publication bias test
abundance_red_all_in_bias_model<-rma.mv(yi,vi,mods = ~-1+sqrt_inv_n_tilda+
                                                        year.c+perc_annual_dist*Functional_group_size+
                                                        I(perc_annual_dist^2)*Functional_group_size, 
                                                        random=list(~1|Site_ID/Study_ID,~1|obsID),
                                                        data=abundance_filtered,
                                                        method="REML",
                                                        test="t")
#suggests that there is a marginal small-study effect, i.e. that effect sizes with larger uncertainty tend to be larger

###############################################################################
#2 - models of heterogeneity in response of diversity of soil and litter fauna#
###############################################################################


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


