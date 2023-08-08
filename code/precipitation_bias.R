#script to examine biases in precipitation

rm(list = ls())


#load the packages
pacman::p_load(raster,tidyverse,ggpubr,geodata,ggbeeswarm,nlme,lme4,bernr,cowplot)

#import site data
site_data<-read.csv("data/fauna_spatial_data.csv")

#summarise data to give numbers of comparisons from each site for each unique 
site_data_unique<-site_data%>%mutate(disturbance_type2=if_else(disturbance_type=="drought",
                                     "Precipitation\nreduction","Precipitation\nincrease"))%>%
  group_by(Study_ID,Site_ID,lat,lon,perc_annual_dist,disturbance_type2)%>%
  summarise(count=length(disturbance_type2))%>%
  filter(!is.na(perc_annual_dist)) #remove rows with no estimate of precipitation change during studies


#get future data and format
#mid-range 2050 scenario
precip_2050<-cmip6_world(model="HadGEM3-GC31-LL",ssp=245,time="2041-2060",var="prec",res=2.5,path="data/spatial_data/climate")
precip_2050_brick<-brick(precip_2050)
precip_2050_total<-calc(precip_2050_brick,fun=sum,filename="data/spatial_data/climate/precip_2050.tif") #sum precipitation for all months

#worst-case  2070 scenario
precip_2070<-cmip6_world(model = "HadGEM3-GC31-LL",ssp = 585,time = "2061-2080",var="prec",res=2.5,path="data/spatial_data/climate")
precip_2070_brick<-brick(precip_2070)
precip_2070_total<-calc(precip_2070_brick,fun=sum,filename="data/spatial_data/climate/precip_2070_worst.tif") #sum precipitation for all months

#get current data and format
precip_present<-worldclim_global("prec",res=2.5,path="data/spatial_data/climate")
precip_present_brick<-brick(precip_present)
precip_present_total<-calc(precip_present_brick,sum,filename="data/spatial_data/climate/precip_current.tif") #sum precipitation for all months

#import saved data
precip_2050_total<-raster("data/spatial_data/climate/precip_2050.tif")
precip_2070_total<-raster("data/spatial_data/climate/precip_2070_worst.tif")
precip_present_total<-raster("data/spatial_data/climate/precip_current.tif")

# extract precip at each site 
coords<-data.frame(lon=site_data_unique$lon, lat=site_data_unique$lat)
coordinates(coords)<-c("lon","lat")

#extract data from raster and append to df
precip_present_extracted <-raster::extract(x=precip_present_total, y=coords)
precip_50_extracted <-raster::extract(x=precip_2050_total, y=coords)
precip_70_extracted <-raster::extract(x=precip_2070_total, y=coords)

#append these to the site dataframe
site_data_unique$precip_present<-precip_present_extracted
site_data_unique$precip_50<-precip_50_extracted
site_data_unique$precip_70<-precip_70_extracted

#calculate precipitation manipulation in mm for studies
site_data_unique$precip_change_study<-site_data_unique$precip_present+(site_data_unique$precip_present*(site_data_unique$perc_annual_dist/100))

#where value is 0 add a 0.1
site_data_unique$precip_change_study_edited<-ifelse(site_data_unique$precip_change_study==0,0.1,site_data_unique$precip_change_study)

#calculate difference in projected vs experimental changes in log response ratios
site_data_unique$precip_exp_proj50<-log(site_data_unique$precip_change_study_edited)-log(site_data_unique$precip_50)
site_data_unique$precip_exp_proj70<-log(site_data_unique$precip_change_study_edited)-log(site_data_unique$precip_70)


#test the difference of this to zero
bias_model<-lme(precip_exp_proj50~disturbance_type2-1,random=~1|Study_ID,data=site_data_unique)
summary(bias_model)

#create new data for prediction
new_data<-data.frame(disturbance_type2=unique(site_data_unique$disturbance_type2))

#predict response
bias_preds<-bolker_ci(bias_model,new_data,pred_int = TRUE,conf_level = 0.95)

#violin plot of this data
ggplot(site_data_unique,aes(x=disturbance_type2,y=precip_exp_proj,fill=disturbance_type2))+
  geom_violin(trim = TRUE)+
  geom_beeswarm(alpha=0.6)+
  geom_hline(yintercept = 0,lty=2)+
  theme_cowplot()+
  scale_fill_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  xlab("Disturbance type")+
  ylab("Difference between study\nand projected precipitation (lnRR)")+
  theme(legend.position = "none") 

ggsave("figures/for_paper/precipitation_bias.png",width = 12,height=10,units = "cm",dpi = 300)

#histogram plot of the data
ggplot(site_data_unique,aes(x=precip_exp_proj,fill=disturbance_type2))+
  geom_histogram()+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  scale_fill_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  xlab("Disturbance type")+
  ylab("Difference between study\nand projected precipitation (lnRR)")+
  theme(legend.position = "none") 

#error plot of this
ggplot(bias_preds,aes(pred,disturbance_type2,colour=disturbance_type2))+
  geom_errorbarh(aes(xmin=ci_l,xmax=ci_h),height=0.2,size=1,alpha=0.8)+
  geom_point(size=4)+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  labs(x="Difference between study and\nprojected precipitation (lnRR)",
       y="Disturbance type")+
  theme(text=element_text(size=12),
      axis.text=element_text(size=10),
      legend.position = "none")

ggsave("figures/for_paper/precipitation_bias_error_plot.png",width = 12,height=10,units = "cm",dpi = 300)



#need to think about how to display all of this best


###########################
#below this is Leo's code##
###########################

###########################################
#1 - get data and format###################
###########################################
?worldclim_global


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
coords<-data.frame(lon=site_data_unique$lon, lat=site_data_unique$lat)
coordinates(coords)<-c("lon","lat")

#extract data from raster and append to df
precip_50 <-data.frame(raster::extract(x=precip_2050, y=coords))
precip_70 <-data.frame(raster::extract(x=precip_2070, y=coords))

col<- precip_50[[1]]
site_data_unique$precip_50<-col
col<-precip_70[[1]]
site_data_unique$precip_70<-col

#convert to a percentage from -100 to 100
site_data$perc_inc <- (site_data$precip_50/site_data$precip)*100
site_data$perc_inc70 <- (site_data$precip_70/site_data$precip)*100
site_data$perc_annual_dist <- ifelse(site_data$perc_annual_dist>100, site_data$perc_annual_dist-100, site_data$perc_annual_dist*(-1))

#create df with just one point for each magnitude of disturbance at each site
Site_unique<-distinct(site_data, Site_ID,perc_annual_dist, .keep_all = TRUE)
fit<-lm(perc_inc70~perc_annual_dist, data = site_data)
summary(fit)

#plot data
a <- ggplot(Site_unique, aes(x = perc_inc, y= perc_annual_dist))+
  geom_point(alpha = 0.2,size =3)+
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
  geom_point(alpha=0.2,size =3)+
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
                bottom = text_grob("Precipitation change projected (%)"))
