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


ALL_df <- read_excel("data/study_data_for_analysis_07_06.xlsx")

#---------------------------------------------------
#adapt dataset for analysis
#---------------------------------------------------

#remove anomalies - outlying disturbance strengths 
ALL_df$perc_annual_dist <- ifelse(ALL_df$perc_annual_dist<500, ALL_df$perc_annual_dist,  NA)
ALL_df$perc_during_dist <- ifelse(ALL_df$perc_during_dist<500, ALL_df$perc_during_dist,  NA)

#change variances from 0 to NA
ALL_df$control_var <- ifelse(ALL_df$control_var !=0, ALL_df$control_var,  NA)
ALL_df$dist_var <- ifelse(ALL_df$dist_var !=0, ALL_df$dist_var, NA)

# change percentage to be continuous from decrease to increase 
ALL_df$perc_annual_dist <- ifelse(ALL_df$perc_annual_dist>100, ALL_df$perc_annual_dist-100, ALL_df$perc_annual_dist*(-1))
ALL_df$perc_during_dist <- ifelse(ALL_df$perc_during_dist>100, ALL_df$perc_during_dist-100, ALL_df$perc_during_dist*(-1))

#convert SE to SD
ALL_df$control_SD<-ifelse(ALL_df$var_type=="SE", ALL_df$control_var*sqrt(ALL_df$control_n), ALL_df$control_var)
ALL_df$dist_SD<-ifelse(ALL_df$var_type=="SE",ALL_df$dist_var*sqrt(ALL_df$dist_n), ALL_df$dist_var)

#---------------------------------------------------
#1. data exploration and effect sizes
#---------------------------------------------------

#log response ratio
ALL_df_rr <- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                    sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = ALL_df)

#impute missing variances with median variance 
ALL_df_rr$vi<-ifelse(is.na(ALL_df_rr$vi),median(ALL_df_rr$vi, na.rm = TRUE), ALL_df_rr $vi)

#ensure complete data set with moderators of interest
#THIS GETS RID OF A LOT OF DATA - SHOULD WE DO THIS OR NOT - IT IS WHAT WE USE TO CARRY OUT ANALYSIS 
ALL_df_rr <- ALL_df_rr[complete.cases(ALL_df_rr[ , c("Functional_group_size","perc_during_dist","perc_annual_dist","time_after_dist_start", "yi")]),]

#subset ALL_df_rr dataset to set variables of interest
#all abundance 
all_ab_rr <- ALL_df_rr %>%filter(broad_outcome == 'abundance')
#all diversity 
all_div_rr <- ALL_df_rr %>%filter(broad_outcome == 'alpha diversity')
#abundance and drought
all_d_ab_rr <- ALL_df_rr %>%filter(broad_outcome == 'abundance')%>% filter(disturbance_type == 'drought') 
#diversity and drought
all_d_div_rr <- ALL_df_rr %>%filter(broad_outcome == 'alpha diversity')%>%filter(disturbance_type == 'drought')
#abundance and precip increases
all_p_ab_rr <- ALL_df_rr %>%filter(broad_outcome == 'abundance')%>%filter(disturbance_type == 'precip_inc') 
#diversity and precip increases
all_p_div_rr <- ALL_df_rr %>%filter(broad_outcome == 'alpha diversity')%>%filter(disturbance_type == 'precip_inc')


#remove cook outliers 
all_ab_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=all_ab_rr)
cooks_a<- cooks.distance(all_ab_m0)
all_ab <- all_ab_rr %>%cbind(cooks_a) %>%filter(cooks_a < 3.0*mean(cooks_a))

all_div_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=all_div_rr)
cooks_b<- cooks.distance(all_div_m0)
all_div <- all_div_rr %>%cbind(cooks_b) %>%filter(cooks_b < 3.0*mean(cooks_b))


#---------------------------------------------------
#1. Plot figures 
#---------------------------------------------------

#Figure 1
all_ab_o1 <-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, random=~1|Site_ID/Study_ID,data=all_ab)
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



