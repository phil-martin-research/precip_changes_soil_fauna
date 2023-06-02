# ORCHARD PLOT FIGURES

rm(list = ls())

library(tidyverse)
library(cowplot)
library(metafor)
library(orchaRd)
library(ggbeeswarm)
library(tidyr)
library(tidymodels)

######################################
#function to calculate I2 for multilevel models

I2_multi<-function(model){
  W <- diag(1/model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
}

#####################################
#notes:
# - I need to add analysis investigating the potential impacts of publication biases
# -I need to add analysis examining the impact of changes in precipitation magnitude


#read in .csv files with soil fauna data
abundance_red<- read_csv("data/abundance_red_data.csv")
abundance_inc<- read_csv("data/abundance_inc_data.csv")
diversity_red<- read_csv("data/diversity_red_data.csv")
diversity_inc<- read_csv("data/diversity_inc_data.csv")

#---------------------------------------------------
#1. data analysis
#---------------------------------------------------

#################################
#analysis of change in abundance#
#################################

########################################
#following reductions in precipitation##
########################################

fauna_ab_red_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_red)#null model

#calculate cook distances
cooks_ab_red_0<-cooks.distance(fauna_ab_red_m0)

#filter out highly influential comparisons
abundance_red_filtered<- abundance_red %>%cbind(cooks_ab_red_0) %>%filter(cooks_ab_red_0 < 3.0*mean(cooks_ab_red_0,na.rm=TRUE))
#this removed 13 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_ab_red_m0_filter<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_red_filtered)#null model
#this shows a 36% reduction with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
abundance_red_appraisal<-abundance_red%>%
                         filter(Validity!="Low validity")

fauna_ab_red_m0_no_low<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_red_appraisal)#model with no low validity studies
#studies with higher robustness tended to show greater reductions in abundance,
#however studies with higher validity had higher reductions in precipitation

#put all this information into a table
#info to include - estimate, se, p val, Q result, I squared

#create loop to do this
model_type<-c("Null model","Outliers removed","Low validity removed")
model_list<-list(fauna_ab_red_m0,fauna_ab_red_m0_filter,fauna_ab_red_m0_no_low)
sensitivity_summary<-data.frame()
for (i in 1:3){
  params<-broom::tidy(model_list[[i]])
  qe<-model_list[[i]]$QE
  qe_p<-model_list[[i]]$QEp
  I2<-I2_multi(model_list[[i]])
  sens_temp<-data.frame(model_type=model_type[i],params,qe,qe_p,I2)
  sensitivity_summary<-rbind(sensitivity_summary,sens_temp)
}


########################################
#following increases in precipitation##
########################################

fauna_ab_inc_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_inc)#null model

#calculate cook distances
cooks_ab_inc_0<-cooks.distance(fauna_ab_inc_m0)

#filter out highly influential comparisons
abundance_inc_filtered<- abundance_inc %>%cbind(cooks_ab_inc_0) %>%filter(cooks_ab_inc_0 < 3.0*mean(cooks_ab_inc_0,na.rm=TRUE))
#this removed 12 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_ab_inc_m0_filter<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_inc_filtered)#null model
#this shows a marginally significant increase of 17%

#run sensitivity analysis based on critical appraisal quality
abundance_inc_appraisal<-abundance_inc%>%
  filter(Validity!="Low validity")

fauna_ab_inc_m0_no_low<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=abundance_inc_appraisal)#model with no low validity studies
#just looking at more robust studies results in a decrease in effect size and loss of significance
#why does this happen?

#put all this information into a table
#info to include - estimate, se, p val, Q result, I squared

#create loop to do this
model_type<-c("Null model","Outliers removed","Low validity removed")
model_list<-list(fauna_ab_inc_m0,fauna_ab_inc_m0_filter,fauna_ab_inc_m0_no_low)
sensitivity_summary_ab_inc<-data.frame()
for (i in 1:3){
  params<-broom::tidy(model_list[[i]])
  qe<-model_list[[i]]$QE
  qe_p<-model_list[[i]]$QEp
  I2<-I2_multi(model_list[[i]])
  sens_temp<-data.frame(model_type=model_type[i],params,qe,qe_p,I2)
  sensitivity_summary_ab_inc<-rbind(sensitivity_summary_ab_inc,sens_temp)
}

#finish this off later!!



##########################################
#analysis of changes in alpha diversity###
##########################################

########################################
#following reductions in precipitation##
########################################

fauna_div_red_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_red)#null model

#calculate cook distances
cooks_div_red_0<-cooks.distance(fauna_div_red_m0)

#filter out highly influential comparisons
diversity_red_filtered<- diversity_red %>%cbind(cooks_div_red_0) %>%filter(cooks_div_red_0 < 3.0*mean(cooks_div_red_0,na.rm=TRUE))
#this removed 3 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_div_red_m0_filter<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_red_filtered)#null model
#this shows an 8% reduction of diversity with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
diversity_red_appraisal<-diversity_red%>%
  filter(Validity!="Low validity")

fauna_div_red_m0_no_low<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_red_appraisal)#model with no low validity studies
#studies with higher robustness tended to show greater reductions in diversity
ggplot(diversity_red,aes(x=Validity,y=perc_annual_dist))+
  geom_violin()
#however studies with higher validity had higher reductions in precipitation

#put all this information into a table
#info to include - estimate, se, p val, Q result, I squared

#create loop to do this
model_type<-c("Null model","Outliers removed","Low validity removed")
model_list<-list(fauna_div_red_m0,fauna_div_red_m0_filter,fauna_div_red_m0_no_low)
div_red_sensitivity_summary<-data.frame()
for (i in 1:3){
  params<-broom::tidy(model_list[[i]])
  qe<-model_list[[i]]$QE
  qe_p<-model_list[[i]]$QEp
  I2<-I2_multi(model_list[[i]])
  sens_temp<-data.frame(model_type=model_type[i],params,qe,qe_p,I2)
  div_red_sensitivity_summary<-rbind(div_red_sensitivity_summary,sens_temp)
}

########################################
#following increases in precipitation##
########################################

fauna_div_inc_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_inc)#null model

#calculate cook distances
cooks_div_inc_0<-cooks.distance(fauna_div_inc_m0)

#filter out highly influential comparisons
diversity_inc_filtered<- diversity_inc %>%cbind(cooks_div_inc_0) %>%filter(cooks_div_inc_0 < 3.0*mean(cooks_div_inc_0,na.rm=TRUE))
#this removed 3 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_div_inc_m0_filter<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_inc_filtered)#null model
#this shows a 6% increase in diversity with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
diversity_inc_appraisal<-diversity_inc%>%
  filter(Validity!="Low validity")

fauna_div_inc_m0_no_low<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=diversity_inc_appraisal)#model with no low validity studies
#not much change in estimate here

#put all this information into a table
#info to include - estimate, se, p val, Q result, I squared

#create loop to do this
model_type<-c("Null model","Outliers removed","Low validity removed")
model_list<-list(fauna_div_inc_m0,fauna_div_inc_m0_filter,fauna_div_inc_m0_no_low)
div_inc_sensitivity_summary<-data.frame()
for (i in 1:3){
  params<-broom::tidy(model_list[[i]])
  qe<-model_list[[i]]$QE
  qe_p<-model_list[[i]]$QEp
  I2<-I2_multi(model_list[[i]])
  sens_temp<-data.frame(model_type=model_type[i],params,qe,qe_p,I2)
  div_inc_sensitivity_summary<-rbind(div_inc_sensitivity_summary,sens_temp)
}



#---------------------------------------------------
#3. Plot figures 
#---------------------------------------------------

#change in abundance

#make my own version of the plots
#bring together predictions from the different models
abun_preds<-distinct(data.frame(predict(fauna_ab_m1_filtered)))
div_preds<-distinct(data.frame(predict(fauna_div_m1_filtered)))
comb_preds<-rbind(abun_preds,div_preds)
comb_preds_2<-data.frame(disturbance=rep(c("Precipitation\nreduction","Precipitation\nincrease"),2),
           outcome=rep(c("Abundance","Alpha diversity"),each=2),
           comb_preds)

#turn predictions into percentages
comb_preds_2<-comb_preds_2%>%
  mutate(perc_pred=(exp(pred)-1)*100,
         per_ci.lb=(exp(ci.lb)-1)*100,
         per_ci.ub=(exp(ci.ub)-1)*100,
         per_pi.lb=(exp(pi.lb)-1)*100,
         per_pi.ub=(exp(pi.ub)-1)*100)

#organise data into one dataset
abundance_filtered2<-abundance_filtered%>%
  select(-cooks_ab_1)
diversity_filtered2<-diversity_filtered%>%
  select(-cooks_div_1)
fauna_filtered<-rbind(abundance_filtered2,diversity_filtered2)

#relabel disturbance types
fauna_filtered<-fauna_filtered%>%
  mutate(disturbance=if_else(disturbance_type=="drought","Precipitation\nreduction","Precipitation\nincrease"),
         outcome=if_else(broad_outcome=="abundance","Abundance","Alpha diversity"))

#make my own orchard plot
ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=fauna_filtered,aes(x=yi,y=disturbance,colour=outcome,group=outcome,size=1/vi),
                   dodge.width = 1,alpha=0.5)+
  geom_errorbarh(data=comb_preds_2,aes(y=disturbance,xmin=pi.lb,xmax=pi.ub,group=outcome),
                 position=position_dodge(width=1),size=1.5,height=0,colour="black",alpha=0.5)+
  geom_errorbarh(data=comb_preds_2,aes(xmin=ci.lb,xmax=ci.ub,y=disturbance,group=outcome),
                 position=position_dodge(width=1),size=3,height=0,colour="black",alpha=0.5)+
  geom_point(data=comb_preds_2,aes(x=pred,y=disturbance,colour=outcome,fill=outcome),
             position=position_dodge(width=1),size=6,shape=21,colour="black")+
  theme_cowplot()+
  scale_fill_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_color_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_size_continuous(range = c(1,10))+
  labs(y="Disturbance type",x="Soil and litter fauna relative to\nno disturbance (log response ratio)")+
  guides(size = "none")+
  theme(text=element_text(size=12),
        axis.text=element_text(size=10))
ggsave("figures/for_paper/abun_div_summary.png",width = 20,height = 14,units = "cm",dpi = 300)
  

#alternative facetted version of the figure
ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=fauna_filtered,aes(x=yi,y=disturbance,colour=disturbance,group=outcome,size=1/vi),
                   dodge.width = 1,alpha=0.5)+
  geom_errorbarh(data=comb_preds_2,aes(y=disturbance,xmin=pi.lb,xmax=pi.ub),
                 position=position_dodge(width=1),size=1.5,height=0,colour="black",alpha=0.8)+
  geom_errorbarh(data=comb_preds_2,aes(xmin=ci.lb,xmax=ci.ub,y=disturbance),
                 position=position_dodge(width=1),size=3,height=0,colour="black",alpha=0.8)+
  geom_point(data=comb_preds_2,aes(x=pred,y=disturbance,colour=disturbance,fill=disturbance),
             position=position_dodge(width=1),size=4,shape=21,colour="black")+
  theme_cowplot()+
  facet_wrap(~outcome,scales = "free_x")+
  scale_fill_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_color_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_size_continuous(range = c(1,10))+
  labs(y="Disturbance type",x="Soil and litter fauna relative to\nno disturbance (log response ratio)")+
  guides(size = "none")+
  theme(text=element_text(size=12),
        axis.text=element_text(size=10),
        legend.position = "bottom",
        legend.justification = "centre")
ggsave("figures/for_paper/abun_div_summary_facet.png",width = 20,height = 14,units = "cm",dpi = 300)



#alternative version using the orchaRd package
#abundance model
orchard_abun_plot<-orchard_plot(fauna_ab_m1_filtered,  group ='Site_ID',mod = 'disturbance_type',
                  data = abundance_filtered, 
                  xlab = "Change in abundance relative\nto baseline (log response ratio)",
                  k = TRUE, g = TRUE,k.pos="left")+
  theme_cowplot()+
  scale_fill_manual(values = c("#fde725","#1f9e89"))+
  scale_color_manual(values = c("#fde725","#1f9e89"))+
  scale_x_discrete(labels=c("Precipitation\ndecrease", "Precipitation\n increase"))+
  annotate("text", x = 2.3, y = 2, label = "+48%")+
  annotate("text", x = 1.3, y = 2, label = "-48%")+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text = element_text(size=10))

#diversity model
orchard_div_plot<-orchard_plot(fauna_div_m1_filtered,  group ='Site_ID',mod = 'disturbance_type',
                                data = diversity_filtered, 
                                xlab = "Change in alpha diversity relative\nto baseline (log response ratio)",
                                k = TRUE, g = TRUE,k.pos="left")+
  theme_cowplot()+
  scale_fill_manual(values = c("#fde725","#1f9e89"))+
  scale_color_manual(values = c("#fde725","#1f9e89"))+
  scale_x_discrete(labels=c("Precipitation\ndecrease", "Precipitation\n increase"))+
  annotate("text", x = 2.3, y = 0.6, label = "+12%")+
  annotate("text", x = 1.3, y = 0.6, label = "-9%")+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text = element_text(size=10))

#combine into one figure
combined_orchard_plots<-plot_grid(orchard_abun_plot,orchard_div_plot,labels = c("(a)","(b)"))
save_plot("figures/for_paper/combined_orchard_plots.png",combined_orchard_plots,base_height = 12,base_width = 20,units="cm")



#######################################################
#Leo's analysis########################################
#######################################################




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



