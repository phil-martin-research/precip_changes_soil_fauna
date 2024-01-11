# script to explore heterogeneity in the effects of precipitation changes on soil and litter fauna 

rm(list = ls())

pacman::p_load(tidyverse,metafor,cowplot,orchaRd,ggbeeswarm,tidyr,ggthemes,sp,broom,lemon,MuMIn,glmulti,PerformanceAnalytics,GGally,gt,geodata,
               ggmap,mapproj,glmulti)


#load data
soil_fauna_df<- read_csv("data/study_data_for_analysis_07_06.csv")
abundance<- read_csv("data/abundance_data.csv")
richness<- read_csv("data/richness_data.csv")
shannon<- read_csv("data/shannon_data.csv")

###############################################################################
#1 - SUMMARY META ANALYSIS.R#
###############################################################################
#this script is produces summary meta-analyses and associated sensitivity tests
#for changes in abundance, taxonomic richness, and shannon diversity of soil and litter fauna
#as a result of increases and decreases in precipitation
######################################
#function to calculate I2 for multilevel models
I2_multi<-function(model){
  W <- diag(1/model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
}

#####################################

#read in .csv files with soil fauna data
abundance<-read.csv("data/abundance_data.csv")
richness<-read.csv("data/richness_data.csv")
shannon<-read.csv("data/shannon_data.csv")
abundance_red<- read.csv("data/abundance_red_data.csv")
abundance_inc<- read.csv("data/abundance_inc_data.csv")
richness_red<- read.csv("data/richness_red_data.csv")
richness_inc<- read.csv("data/richness_inc_data.csv")
shannnon_red<- read.csv("data/shannon_red_data.csv")
shannnon_inc<- read.csv("data/shannon_inc_data.csv")

######################################################
#1. meta-analyses#####################################
######################################################

fauna_list<-list(abundance_red,abundance_inc,richness_red,richness_inc,shannnon_red,shannnon_inc)
outcomes<-c("Abundance","Abundance","Taxonomic richness","Taxonomic richness","Shannon-Wiener index","Shannon-Wiener index")
disturbances <- rep(c("Precipitation\nreduction", "Precipitation\nincrease"), times = 3)
sensitivity_summary<-data.frame()
prediction_summary<-data.frame()
#loop through all the different stages for the meta-analyses
for (i in 1:length(fauna_list)){
  #subset list to give relevant dataframe
  temp_df<-fauna_list[[i]]
  #run null model
  m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=temp_df)
  #remove comparisons that fail geary's test
  no_geary<-temp_df%>%
    mutate(geary_test=if_else(is.na(geary_test),"Not needed",geary_test))%>%
    filter(geary_test!="fail")
  #null model excluding the studies that fail Geary's test
  m0_no_geary<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=no_geary)
  #remove studies that have low validity as assessed by critical appraisal
  temp_df_appraisal<-temp_df%>%
    filter(Validity!="Low validity")
  #run model with no low validity studies
  m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=temp_df_appraisal)
  #put all this information into a table
  #info to include - estimate, se, p val, Q result, I squared
  #create loop to do this
  model_type<-c("Null model","Failed Geary test","Low validity removed")
  model_list<-list(m0,m0_no_geary,m0_no_low)
  for (y in 1:3){
    params<-broom::tidy(model_list[[y]])
    qe<-model_list[[y]]$QE
    qe_p<-model_list[[y]]$QEp
    I2<-I2_multi(model_list[[y]])
    k<-model_list[[y]]$k.all
    sens_temp<-data.frame(disturbance=disturbances[i],outcome=outcomes[i],
                          model_type=model_type[y],params,k,qe,qe_p,I2)
    sensitivity_summary<-rbind(sensitivity_summary,sens_temp)
  }
  #create predictions based on models
  temp_pred<-distinct(data.frame(predict(m0)))
  temp_pred2<-data.frame(disturbance=disturbances[i],detailed_outcome=outcomes[i],temp_pred)
  prediction_summary<-rbind(prediction_summary,temp_pred2)
}

#combine tables of sensitivity analyses into one big table
sensitivity_table <- sensitivity_summary%>%
  mutate(across(where(is.numeric), round, 3))%>%
  mutate(I2=round(I2,0))%>%
  gt()

#export this to a word file
sensitivity_table%>%gtsave("figures/for_paper/summary_sensitivity_table.docx")


#convert predictions to percentages
prediction_summary_perc<-prediction_summary%>%
  mutate(perc_pred=(exp(pred)-1)*100,
         per_ci.lb=(exp(ci.lb)-1)*100,
         per_ci.ub=(exp(ci.ub)-1)*100,
         per_pi.lb=(exp(pi.lb)-1)*100,
         per_pi.ub=(exp(pi.ub)-1)*100)

#change names of different disturbances and outcomes

#---------------------------------------------------
#3. Plot figures 
#---------------------------------------------------

#produce my own version of an orchard plot

#first produce a plot for abundance

#Organise data into one dataset

fauna_data<-rbind(abundance,richness,shannon)

#relabel disturbance types
fauna_data<-fauna_data%>%
  mutate(disturbance=if_else(disturbance_type=="drought","Precipitation\nreduction","Precipitation\nincrease"),
         detailed_outcome=if_else(detailed_outcome=="abundance","Abundance",
                                  if_else(detailed_outcome=="taxonomic richness","Taxonomic richness","Shannon-Wiener index")))%>%
  filter(lnrr_laj>-3)

#facetted version of the figure
facet_summary_plot<-ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=subset(fauna_data,lnrr_laj<2.5),
                   aes(x=lnrr_laj,y=disturbance,colour=disturbance,group=detailed_outcome,size=1/v_lnrr_laj),
                   dodge.width = 0.01,alpha=0.5)+
  geom_errorbarh(data=prediction_summary_perc,aes(y=disturbance,xmin=pi.lb,xmax=pi.ub),
                 position=position_dodge(width=1),linewidth=1,height=0,colour="black",alpha=0.8)+
  geom_errorbarh(data=prediction_summary_perc,aes(xmin=ci.lb,xmax=ci.ub,y=disturbance),
                 position=position_dodge(width=1),linewidth=2,height=0,colour="black",alpha=0.8)+
  geom_point(data=prediction_summary_perc,aes(x=pred,y=disturbance,colour=disturbance,fill=disturbance),
             position=position_dodge(width=1),size=4,shape=21,colour="black")+
  theme_cowplot()+
  facet_rep_wrap(~detailed_outcome,scales = "free_x",repeat.tick.labels = FALSE,ncol=3)+
  scale_fill_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_size_continuous(range = c(1,10))+
  labs(y="Disturbance type",x="Change in soil and litter fauna outcome (log response ratio)")+
  guides(size = "none")+
  theme(text=element_text(size=20),
        axis.text=element_text(size=14),
        legend.position = "none",
        legend.justification = "centre",
        strip.background =element_rect(fill="grey",color = "grey", linewidth = 1),
        strip.text = element_text(face = "bold",margin = unit(rep(2, 4), "pt")),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
facet_summary_plot






#create data for the labelling using number of studies and number of comparisons
studies_label<-fauna_data%>%
  group_by(disturbance,detailed_outcome)%>%
  summarise(k_count=length(lnrr_laj),
            studies=n_distinct(Study_ID),
            studies_annotation=paste("k = ",k_count," (",studies,")",sep=""))
studies_label$lnrr_laj<-c(-1.5,-1.5,-1.5,-1.5,-1.5,-1.5)

#create data for labelling using effect sizes in percentages etc
effect_size_label<-prediction_summary_perc%>%
  mutate(perc_change=round(perc_pred,0),
         change_label=if_else(perc_change<0,
                              paste(perc_change,"%",sep=""),
                              paste("+",perc_change,"%",sep="")),
         change_label=if_else(sign(per_ci.lb)==sign(per_ci.ub),paste(change_label,"*",sep = ""),change_label))
effect_size_label$lnrr_laj<-c(1.3,1.3,0.1,0.1,-.5,-.5)

#add these to the plot
facet_summary_plot_with_label1<-facet_summary_plot+
  geom_text(data=studies_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=studies_annotation),
            hjust   = 0.5,
            vjust   = -2.2,
            size = 5)


facet_summary_plot_with_label2<-facet_summary_plot_with_label1+
  geom_text(data=effect_size_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=change_label),
            hjust   = 0,
            vjust   = -2.2,
            size=5)
facet_summary_plot_with_label2


#add labels for figures
facet_summary_plot_with_label3<-egg::tag_facet(facet_summary_plot_with_label2)+
  theme(strip.text = element_text(face = "bold",margin = unit(rep(2, 4), "pt")),
        strip.background = element_rect(fill="white",color = "black", size = 1))
facet_summary_plot_with_label3

ggsave("figures/for_poster/abun_div_summary_facet.png",facet_summary_plot_with_label3,width = 13,height = 20,units = "cm",dpi = 300)
ggsave("figures/for_poster/abun_div_summary.png",facet_summary_plot_with_label2,width = 30,height = 10,units = "cm",dpi = 300)






###############################################################################
#2 - HETEROGENEITY ANALYSIS#
###############################################################################
#------------------------------------------------------------------------------
#1 - models of heterogeneity in response of abundance of soil and litter fauna#
#------------------------------------------------------------------------------
#best model is M12
M12<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size.y+year.c+sqrt_inv_n_tilda-1,random=~1|Study_ID/Site_ID/obsID,data=abundance_complete)

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
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 15),  # Adjust axis title text size
        strip.text = element_text(size = 17), 
        axis.text = element_text(size = 12))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)
#+  geom_text(data=sample_size_label,aes(x=150,y=4,label=k_label),colour="black")

#save plot
ggsave("figures/for_poster/abundance_precip_size.png",width = 20,height = 8,units = "cm",dpi = 300)

###############################################################################
#2 - models of heterogeneity in response of diversity of soil and litter fauna#
###############################################################################

#2.1 - taxonomic diversity

#save best model formula
rich_M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist,random=~1|Study_ID/Site_ID/obsID,data=richness_complete)
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
ggsave("figures/for_poster/richness_precip.png",width = 12,height = 8,units = "cm",dpi = 300)



#2.2 - shannon diversity
# results offer weak support for the hypothesis that there has been a change in impact of drought studies over time

#plot figure showing these results

#save model formula
shannon_M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+year.c,random=~1|Study_ID/Site_ID/obsID,data=shannon_complete)
M5_formula<-(~perc_annual_dist+year.c)


#and now a figure for precipitation impact
#create new dataset for predictions for year
new_data_shannon_precip<-data.frame(perc_annual_dist=seq(min(shannon_complete$perc_annual_dist),max(shannon_complete$perc_annual_dist),1),
                                    year.c=mean(shannon_complete$year.c))

#create a model matrix and remove the intercept
predgrid_shannon_precip<-model.matrix(M5_formula,data=new_data_shannon_precip)[,-1]

#predict onto the new model matrix
mypreds_shannon_precip<-data.frame(predict.rma(shannon_M5,newmods=predgrid_shannon_precip))

#attach predictions to variables for plotting
new_data_shannon_precip<-cbind(new_data_shannon_precip, mypreds_shannon_precip[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot the results of the predictions
shannon_precip_figure<-ggplot(new_data_shannon_precip,aes(perc_annual_dist,y=pred))+
  geom_line(colour="#fcba03")+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA,fill="#fcba03")+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA,fill="#fcba03")+
  geom_point(data=shannon_complete,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),colour="#fcba03",alpha=0.2)+
  theme_cowplot()+
  labs(x="Change in precipitation (%)",
       y="Change in soil fauna Shannon\ndiversity(log response ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)
shannon_precip_figure


save_plot("figures/for_poster/shannon_change.png",shannon_precip_figure,base_height = 10,base_width = 18,units="cm",dpi=300)






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
bias_preds$perc_pred<-(exp(bias_preds$pred)-1)*100
bias_preds$perc_pred_label<-paste(round(bias_preds$perc_pred),"%",sep="")

#error plot of this
precip_bias_plot<-ggplot(bias_preds,aes(pred,disturbance_type2,colour=disturbance_type2))+
  geom_errorbarh(aes(xmin=ci_l,xmax=ci_h),height=0.2,linewidth=1,alpha=0.8)+
  geom_point(size=4)+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  labs(x="Difference between study and\nprojected precipitation (lnRR)",
       y="Disturbance type")+
  geom_text(aes(x = pred, y = disturbance_type2, label = perc_pred_label),
            hjust = 0.3, vjust = -2, colour = "black", size = 6) +  # Adjust text size here
  theme(text=element_text(size=14),
        axis.text=element_text(size=12),
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
        text=element_text(size=14),
        axis.text=element_text(size=12))

#combine the two plots
bias_plots<-plot_grid(precip_bias_plot,plot_bias_plot, ncol = 1)
bias_plots
save_plot("figures/for_poster/bias_plots.png",bias_plots,base_height = 17,base_width = 10,units="cm",dpi=300)
save_plot("figures/for_poster/bias_plots.pdf",bias_plots,base_height =17,base_width = 10,units="cm",dpi=300)


