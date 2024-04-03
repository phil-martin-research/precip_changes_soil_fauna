#script to examine biases in precipitation

rm(list = ls())


#load the packages
pacman::p_load(raster,tidyverse,ggpubr,geodata,ggbeeswarm,nlme,lme4,bernr,cowplot,scales,ggtext)

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

#create string for all the different climate models

climate_models<-c("ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", 
                  "BCC-CSM2-MR", "CanESM5", "CanESM5-CanOE", "CMCC-ESM2", 
                  "CNRM-CM6-1", "CNRM-CM6-1-HR", "CNRM-ESM2-1", 
                  "EC-Earth3-Veg", "EC-Earth3-Veg-LR", "FIO-ESM-2-0", 
                  "GFDL-ESM4", "GISS-E2-1-G", "GISS-E2-1-H", "HadGEM3-GC31-LL", 
                  "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC-ES2L", 
                  "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL")





#get current data and format
precip_present<-worldclim_global("prec",res=2.5,path="data/spatial_data/climate")
precip_present_brick<-brick(precip_present)
precip_present_total<-raster::calc(precip_present_brick,sum,filename="data/spatial_data/climate/precip_current.tif") #sum precipitation for all months

# extract current precipitation at each site 
coords<-data.frame(lon=site_data_unique$lon, lat=site_data_unique$lat)
coordinates(coords)<-c("lon","lat")

precip_present_total<-raster("data/spatial_data/climate/precip_current.tif")
precip_present_extracted <-raster::extract(x=precip_present_total, y=coords)
site_data_unique$precip_present<-precip_present_extracted


?cmip6_world

for (i in 1:3){
  precip_2050<-cmip6_world(model=climate_models[[i]],ssp=245,time="2041-2060",var="prec",res=2.5,path=tempdir())
  precip_2050_brick<-brick(precip_2050)
  precip_2050_total<-raster::calc(precip_2050_brick,fun=sum,filename=paste("data/spatial_data/climate/scenario_sums/precip_2050_",
                                                                           climate_models[[i]],".tif",sep = ""),overwrite=TRUE) #sum precipitation for all months
}

plot(precip_2050)

cmip6_world(model=climate_models[[i]],ssp=245,time="2041-2060",var="prec",res=2.5,path="data/spatial_data/climate",lon=116.972222,lat=4.745833)

#run a loop to:
#1 - import scenario data
#2 - extract data from each scenario
#3 - add column to dataframe for each scenario
#4 - add column with model name to each dataframe
##5 - stack dataframes on top of each other


#import saved data

#list of all scenarios
scenario_list<-list.files("data/spatial_data/climate/scenario_sums/",pattern=".tif")

for (i in 1:length(scenario_list)){
  i<-1
  precip_2050_total<-raster(paste("data/spatial_data/climate/scenario_sums/",scenario_list[[i]],sep=""))
  precip_50_extracted <-raster::extract(x=precip_2050_total, y=coords)
  site_data_unique$precip_50<-precip_50_extracted
  #where value is 0 add a 0.1
  site_data_unique$precip_change_study_edited<-ifelse(site_data_unique$perc_annual_dist==0,0.1,site_data_unique$perc_annual_dist)
  #calculate difference in projected vs experimental changes in log response ratios
  site_data_unique$precip_exp_proj50<-log(site_data_unique$precip_change_study_edited)-log(site_data_unique$precip_50)
}

plot(precip_2050_total)

precip_2050_total<-raster("data/spatial_data/climate/scenario_sums/precip_2050_ACCESS-CM2.tif")



precip_2070_total<-raster("data/spatial_data/climate/precip_2070_worst.tif")








#test the difference of this to zero
bias_model<-lme(precip_exp_proj50~disturbance_type2-1,random=~1|Study_ID,data=site_data_unique)
summary(bias_model)

#create new data for prediction
new_data<-data.frame(disturbance_type2=unique(site_data_unique$disturbance_type2))

#predict response
bias_preds<-bolker_ci(bias_model,new_data,pred_int = TRUE,conf_level = 0.95)
bias_preds$perc_pred<-(exp(bias_preds$pred)-1)*100
bias_preds$perc_pred_label<-paste(round(bias_preds$perc_pred),c("%*","%"),sep="")

#error plot of this
precip_bias_plot<-ggplot(bias_preds,aes(pred,disturbance_type2,colour=disturbance_type2))+
  geom_errorbarh(aes(xmin=ci_l,xmax=ci_h),height=0.2,linewidth=1,alpha=0.8)+
  geom_point(size=4)+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  labs(x="Difference between study and\nprojected precipitation (lnRR)",
       y="Disturbance type")+
  geom_text(aes(x=c(-2,2),y=disturbance_type2,label=perc_pred_label),hjust   = 0.3,vjust   = -2,colour="black")+
  theme(text=element_text(size=12),
      axis.text=element_text(size=10),
      legend.position = "none")
precip_bias_plot
ggsave("figures/for_paper/precipitation_bias_error_plot.png",width = 12,height=10,units = "cm",dpi = 300)



#bias associated with plot area and body size
site_data<-read.csv("data/site_data_2023_06_07.csv")
fact_table<-read.csv("data/fact_table_2023_08_06.csv")

site_bias_data<-fact_table%>%left_join(site_data,"Site_ID")%>%
                             filter(use_for_first_analysis==TRUE,!is.na(Functional_group_size))%>%
                             distinct(Site_ID,plot_area.y,perc_annual_dist,Functional_group_size,exp_obs)%>%
                             filter(Functional_group_size!="",exp_obs=="experimental")


site_bias_data%>%
  group_by(Functional_group_size)%>%
  summarise(med_area=median(plot_area.y,na.rm=TRUE))

names(site_bias_data)

plot_bias_plot<-ggplot(site_bias_data,aes(Functional_group_size,plot_area.y,fill=Functional_group_size,colour=Functional_group_size))+
  geom_boxplot(alpha=0.5)+
  theme_cowplot()+
  scale_x_discrete(limits=rev)+
  labs(x="Body size",y=bquote("Plot area "(m^2)))+
  scale_color_manual(values = c("#e0b500","#c3386b","#02475f"))+
  scale_fill_manual(values = c("#e0b500","#c3386b","#02475f"))+
  theme(legend.position="none",
        text=element_text(size=12),
        axis.text=element_text(size=10))

#combine the two plots
bias_plots<-plot_grid(precip_bias_plot,plot_bias_plot,labels = c("(a)","(b)"))
save_plot("figures/for_paper/bias_plots.png",bias_plots,base_height = 10,base_width = 20,units="cm",dpi=300)
save_plot("figures/for_paper/bias_plots.pdf",bias_plots,base_height = 10,base_width = 20,units="cm",dpi=300)
