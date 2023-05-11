# ORCHARD PLOT FIGURES

rm(list = ls())

library(colorblindr)
library(scico)
library(ggpubr)
library(cowplot)
library(readxl)
library(ggplot2)
library(metafor)
library(dplyr)
library(raster)
library(maps)
library(mapdata)
library(RColorBrewer)
library(cowplot)
library(orchaRd)
require(scales)
require(rgdal)
require(readr)
require(dismo)
require(maptools)
require(CoordinateCleaner)
require(ggsn)
require(rgeos)
require(sf)
require(ggsn)
require(ggspatial)
library(readr)


#read in .csv file with soil fauna data
soil_fauna_df<- read_csv("data/study_data_for_analysis_07_06.csv")

#---------------------------------------------------
#format dataset for analysis
#---------------------------------------------------

#clean dataset
soil_fauna_df <- soil_fauna_df %>%
  mutate(
    #remove disturbance strengths that represent >500% increases in precipitation
    perc_annual_dist = if_else(perc_annual_dist < 500, perc_annual_dist, NA_real_),
    perc_during_dist = if_else(perc_during_dist < 500, perc_during_dist, NA_real_),
    #change all varainces that are 0 to NA
    control_var = if_else(control_var !=0, control_var, NA_real_),
    dist_var = if_else(dist_var !=0, dist_var, NA_real_),
    #alter format of precipitation change to represent continuous change from increase to decrease
    perc_annual_dist=if_else(perc_annual_dist>100, perc_annual_dist-100, perc_annual_dist*(-1)),
    perc_during_dist=if_else(perc_during_dist>100, perc_during_dist-100, perc_during_dist*(-1)),
    #convert SE to SD
    control_SD=if_else(var_type=="SE", control_var*sqrt(control_n), control_var),
    dist_SD=if_else(var_type=="SE",dist_var*sqrt(dist_n), dist_var)
  )

#check to see if any of the means are equal to zero
#control group
soil_fauna_df%>%
  group_by(control_av)%>%
  summarise(length(control_av))
#disturbance group
soil_fauna_df%>%
  group_by(disturbance_av)%>%
  summarise(length(disturbance_av))

#there are only 5 data points where control or disturbance mean are equal to 0
#so we will exclude these
soil_fauna_df<-soil_fauna_df%>%
  filter(control_av>0&disturbance_av>0)

#impute missing SD values based on the coefficient of variation
#and missing sample sizes based on median sample sizes
control_cv<-median(soil_fauna_df$control_SD/soil_fauna_df$control_av,na.rm = TRUE)
disturbance_cv<-median(soil_fauna_df$dist_SD/soil_fauna_df$disturbance_av,na.rm = TRUE)
med_control_n<-median(soil_fauna_df$control_n,na.rm = TRUE)
med_dist_n<-median(soil_fauna_df$dist_n,na.rm = TRUE)

soil_fauna_df<-soil_fauna_df%>%
  mutate(control_SD=if_else(is.na(control_SD),control_cv*control_av,control_SD),
         dist_SD=if_else(is.na(dist_SD),disturbance_cv*disturbance_av,dist_SD),
         control_n=if_else(is.na(control_n),med_control_n,control_n),
         dist_n=if_else(is.na(dist_n),med_dist_n,dist_n))

#calculate log response ratio
soil_fauna_rr<- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                    sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = soil_fauna_df)


#---------------------------------------------------
#1. data exploration and effect sizes
#---------------------------------------------------


#subset dataset to set variables of interest
#all abundance 
fauna_ab<- soil_fauna_rr %>%filter(broad_outcome == 'abundance')
#all diversity 
fauina_div<- soil_fauna_rr %>%filter(broad_outcome == 'alpha diversity')

####I don´t think we need the subset below this######

#abundance and drought
all_d_ab_rr <- ALL_df_rr %>%filter(broad_outcome == 'abundance')%>% filter(disturbance_type == 'drought') 
#diversity and drought
all_d_div_rr <- ALL_df_rr %>%filter(broad_outcome == 'alpha diversity')%>%filter(disturbance_type == 'drought')
#abundance and precip increases
all_p_ab_rr <- ALL_df_rr %>%filter(broad_outcome == 'abundance')%>%filter(disturbance_type == 'precip_inc') 
#diversity and precip increases
all_p_div_rr <- ALL_df_rr %>%filter(broad_outcome == 'alpha diversity')%>%filter(disturbance_type == 'precip_inc')


#remove cook outliers 

#analysis of change in abundance
fauna_ab_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=fauna_ab)#null model
fauna_ab_m1 <-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, random=~1|Site_ID/Study_ID,data=fauna_ab)#impact of drought vs increases

#calculate cook distances
cooks_ab_0<-cooks.distance(fauna_ab_m0)
cooks_ab_1<-cooks.distance(fauna_ab_m1)
#combine into one dataframe
c_dists<-data.frame(cooks_ab_0,cooks_ab_1)
#compare cook distances for different models
ggplot(c_dists,aes(cooks_ab_0,cooks_ab_1))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

#there are some cooks distances for model 1 that are much larger than for the null model
#this means that we should filter out high cooks distances for model 1 and then rerun the model

fauna_ab_filtered<- fauna_ab %>%cbind(cooks_ab_1) %>%filter(cooks_ab_1 < 3.0*mean(cooks_ab_1))

#rerun analysis of impact of decreases vs decreases
fauna_ab_m1_filtered<-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, 
                             random=~1|Site_ID/Study_ID,data=fauna_ab_filtered)


fauna_ab_filtered$disturbance_type

all_div_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=all_div_rr)
cooks_b<- cooks.distance(all_div_m0)
all_div <- all_div_rr %>%cbind(cooks_b) %>%filter(cooks_b < 3.0*mean(cooks_b))





#---------------------------------------------------
#1. Plot figures 
#---------------------------------------------------

#Figure 1

all_div_o1 <-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, random=~1|Site_ID/Study_ID,data=all_div)
c <- orchard_plot(all_ab_o1,  group ='Site_ID',mod = 'disturbance_type',
                  data = all_ab, xlab = "Abundance relative to undisturbed soil \nlog(Response ratio",
                  k = FALSE, g = FALSE)+
  ggtitle('Abundance')+
  theme_cowplot()+
  scale_fill_manual(values = c("#fde725","#1f9e89"))+
  scale_color_manual(values = c("#fde725","#1f9e89"))+
  scale_x_discrete(labels=c("Precipitation\ndecrease", "Precipitation\n increase"))

d <- orchard_plot(all_div_o1, group ='Site_ID',mod = 'disturbance_type',
                  data = all_div, xlab = "Diversity relative to undisturbed soil \nlog(Response ratio",
                  k = FALSE, g = FALSE)+
  ggtitle('Diversity')+
  theme_cowplot()+
  scale_fill_manual(values = c("#fde725","#1f9e89"))+
  scale_color_manual(values = c("#fde725","#1f9e89"))+
  scale_x_discrete(labels=c("Precipitation\ndecrease", "Precipitation\n increase"))

ggarrange(c +rremove("legend"), d +rremove("legend")+ rremove("ylab"), nrow = 2)


#Figure 2. Diversity magnitude 
all_div_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=all_div_rr)
cooks_c<- cooks.distance(all_div_m0)
cooks_div_rr <- all_div_rr %>%cbind(cooks_c) %>%filter(cooks_c < 3.0*mean(cooks_c))

div_m1<-rma.mv(yi,vi,mods = ~perc_annual_dist, random=~1|Site_ID/Study_ID,data=cooks_div_rr)
pred_div_m1 <- predict(div_m1, newmods=c(-100:80), addx=TRUE)
pred_div_m1_df<- as.data.frame(pred_div_m1)

# plot cooks_div strength
div_f1 <- ggplot()+
  geom_point(data = cooks_div_rr, aes(x = perc_annual_dist , y = yi, size = (1/vi)), alpha =0.5)+
  scale_size_continuous(range = c(3, 8))+
  geom_ribbon(data = pred_div_m1_df, 
              aes(x = X.perc_annual_dist, y = pred, ymin = ci.lb, ymax = ci.ub), alpha = 0.2)+
  geom_line(data = pred_div_m1_df, 
            aes(x = X.perc_annual_dist, y = pred), size=1)+
  theme_cowplot(font_size = 18)+
  ylab("Diversity relative to undisturbed soil \nlog(response ratio)")+
  xlab("Percentage of mean annual precipitation")
div_f1


#Figure 3. Precipitation biases compared to CMIP5 predictions. 

#get environmental data
env_2050 <- getData('CMIP5', var="bio", res=2.5, model="HE", year=50, rcp=85)
env_2070 <- getData('CMIP5', var="bio", res=2.5, model="HE", year=70, rcp=85)
env <- getData("worldclim", var="bio", res=2.5)

# name data 
bioclim_names <-
  c("Annual_Mean_Temp","Mean_Diurnal_Range","Isothermality","Temp_Seasonality",
    "Max_Temp_Warmest Month", "Min_Temp_Coldest_Month", "Temp_Annual_Range",
    "Mean_Temp_Wettest_Quarter","Mean_Temp_Driest_Quarter","Mean_Temp_Warmest_Quarter",
    "Mean_Temp_Coldest_Quarter","Annual_Precip","Precip_Wettest_Month",
    "Precip_Driest_Month","Precip_Seasonality","Precip_Wettest_Quarter",
    "Precip_Driest_Quarter","Precip_Warmest_Quarter","Precip_Coldest_Quarter"
  )

names(env) <- bioclim_names
names(env_2050)<- bioclim_names

#find difference in annual precipitation 
precip_2050 <- env_2050[[12]]-env[[12]]
precip_2070 <- env_2070[[12]]-env[[12]]
perc_2050 <- (precip_2050/env[[12]])*100

# extract temp and precip at each site 
coords<-data.frame(lon=ALL_df$lon, lat=ALL_df$lat)
coordinates(coords)<-c("lon","lat")

#extract data from raster and append to df
precip_50 <-data.frame(raster::extract(x=precip_2050, y=coords))
precip_70 <-data.frame(raster::extract(x=precip_2070, y=coords))

col<- precip_50[[1]]
ALL_df$precip_50<-col
col<-precip_70[[1]]
ALL_df$precip_70<-col

#convert to a percentage from -100 to 100
ALL_df$perc_inc <- (ALL_df$precip_50/ALL_df$precip)*100
ALL_df$perc_inc70 <- (ALL_df$precip_70/ALL_df$precip)*100
ALL_df$perc_annual_dist <- ifelse(ALL_df$perc_annual_dist>100, ALL_df$perc_annual_dist-100, ALL_df$perc_annual_dist*(-1))

#create df with just one point for each magnitude of disturbance at each site
Site_unique<-distinct(ALL_df, Site_ID,perc_annual_dist, .keep_all = TRUE)
fit<-lm(perc_inc70~perc_annual_dist, data = Site_unique)
summary(fit)

#plot data
a <- ggplot(Site_unique, aes(x = perc_inc, y= perc_annual_dist))+
  geom_point(aes(alpha = 0.2),size =3)+
  xlim(-100,100)+
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  ggtitle("A")+
  xlab("Projected precipitation change (%)")+
  ylab("Experimental precipitation change (%)")+
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = .2)+
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = .2)+
  annotate(geom="text", x=-70, y=130, label="2050",color="black")+
  geom_abline(intercept = 0, slope = 1)+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)

b<- ggplot(Site_unique, aes(x= perc_inc70, y = perc_annual_dist))+
  geom_point(aes(alpha = 0.2),size =3)+
  xlim(-100,100)+
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  ggtitle("B")+
  xlab("Projected precipitation change (%)")+
  ylab("Experimental precipitation change (%)")+
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, 
           alpha = .2)+
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, 
           alpha = .2)+
  annotate(geom="text", x=-70, y=130, label="2070",color="black")+
  geom_abline(intercept = 0, slope = 1)+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)


figure<- ggarrange(a+ rremove("ylab")+ rremove("xlab"),b + rremove("xlab")+ rremove("ylab"),nrow = 1,common.legend = TRUE, legend = "right")
annotate_figure(figure, left = text_grob("Precipitation change experimental (%)", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = text_grob("Precipitation change projected (%)", gp = gpar(cex = 1.3)))



