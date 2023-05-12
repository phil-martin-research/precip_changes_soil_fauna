# ORCHARD PLOT FIGURES

rm(list = ls())

library(tidyverse)
library(cowplot)
library(metafor)
library(orchaRd)
library(ggbeeswarm)
library(tidyr)


#read in .csv file with soil fauna data
soil_fauna_df<- read_csv("data/study_data_for_analysis_07_06.csv")

#---------------------------------------------------
#1 - format dataset for analysis
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

#extract year of study
soil_fauna_df$study_year<-parse_number(soil_fauna_df$Study_ID,trim_ws = TRUE)

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

#subset dataset to get variables of interest
#all abundance 
fauna_ab<- soil_fauna_rr %>%filter(broad_outcome == 'abundance')
#all diversity 
fauna_div<- soil_fauna_rr %>%filter(broad_outcome == 'alpha diversity')

######################################################
#2 - test of publication bias#########################
######################################################

#subset to give only studies of each precipitiation change and outcome combination
fauna_ab_red<-fauna_ab %>%filter(disturbance_type == 'drought')
fauna_ab_inc<-fauna_ab %>%filter(disturbance_type == 'precip_inc')
fauna_div_red<-fauna_div %>%filter(disturbance_type == 'drought')
fauna_div_inc<-fauna_div %>%filter(disturbance_type == 'precip_inc')

#precipitation reduction and abundance
par(mfrow = c(4, 2))
funnel(fauna_ab_red$yi, fauna_ab_red$vi, yaxis="sei",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)", digits = 2, las = 1) 
funnel(fauna_ab_red$yi, fauna_ab_red$vi, yaxis="seinv",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)",  digits = 2, las = 1) 
#there seems to be a slight bias towards more negative effects

#precipitation reduction and diversity
funnel(fauna_ab_inc$yi, fauna_ab_inc$vi, yaxis="sei",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)", digits = 2, las = 1) 
funnel(fauna_ab_inc$yi, fauna_ab_inc$vi, yaxis="seinv",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)",  digits = 2, las = 1) 
#again a slight bias towards more negative effects

#precipitation increases and abundance
funnel(fauna_ab_red$yi, fauna_ab_red$vi, yaxis="sei",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)", digits = 2, las = 1) 
funnel(fauna_ab_red$yi, fauna_ab_red$vi, yaxis="seinv",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)",  digits = 2, las = 1) 
#less bias obvious here

#precipitation increase and diversity
funnel(fauna_div_inc$yi, fauna_div_inc$vi, yaxis="sei",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)", digits = 2, las = 1) 
funnel(fauna_div_inc$yi, fauna_div_inc$vi, yaxis="seinv",
       #xlim = c(-3, 3),
       xlab = "Effect size (lnRR)",  digits = 2, las = 1) 
#little evidence of bias


#---------------------------------------------------
#3. data analysis
#---------------------------------------------------

#################################
#analysis of change in abundance#
#################################

fauna_ab_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=fauna_ab)#null model
fauna_ab_m1 <-rma.mv(yi,vi,mods = ~disturbance_type-1, random=~1|Site_ID/Study_ID,data=fauna_ab)#impact of drought vs increases

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
fauna_ab_filtered<- fauna_ab %>%cbind(cooks_ab_1) %>%filter(cooks_ab_1 < 3.0*mean(cooks_ab_1,na.rm=TRUE))
#this removed 8 comparisons

#rerun analysis of impact of decreases vs decreases
fauna_ab_m1_filtered<-rma.mv(yi,vi,mods = ~disturbance_type-1, 
                             random=~1|Site_ID/Study_ID,data=fauna_ab_filtered)

#this shows a 48% reduction with precipitation decreases and a 49% increase for precipitation increases

###########################################################
#tests of publication bias for abundance###################
###########################################################

#this section tests for publication bias following the recommendations of Nakagawa et al (2021) https://doi.org/10.1111/2041-210X.13724

#1 - precipitation reductions

#create a variable for the standard error of each effect size
fauna_ab_red$sei <- sqrt(fauna_ab_red$vi)

#run a meta-analytical model of impacts of preciptation reduction
fauna_ab_m0_reduction<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=fauna_ab_red)

#produce a funnel plot
par(mfrow=c(1,1))
funnel(fauna_ab_red$yi,
       fauna_ab_red$vi,
       yaxis="seinv",
       ylab="Precision (1/SE)",
       xlab = "Effect size (lnRR)")
#this generally looks fine, although there is some bias towards negative values

#calculate effective sample size 
fauna_ab_red$e_n<-with(fauna_ab_red,(4*(control_n*dist_n)) / (control_n + dist_n))
#plot funnel plot with effective sample size
funnel(fauna_ab_red$yi, fauna_ab_red$vi, ni = fauna_ab_red$e_n, yaxis="ni",
       #xlim = c(-3, 3),
       ylab = "Effective sample size",
       xlab = "Effect size (lnRR)") 



#use multilevel models to examine publication biases

#meta-regression with SE
abundance_red_se_bias_model<-rma.mv(yi,vi,mods = ~1+sei, 
                                random=~1|Site_ID/Study_ID,
                                data=fauna_ab_red,
                                method="REML",
                                test="t")
#this indicates that there is little evidence of small-study effects

#meta-regression with year of publication

# mean-centering year of publication to help with interpretation
fauna_ab_red$year.c <- as.vector(scale(fauna_ab_red$study_year, scale = F))

#meta-regression of effect of publication year
abundance_red_year_bias_model<-rma.mv(yi,vi,mods = ~1+year.c, 
                                random=~1|Site_ID/Study_ID,
                                data=fauna_ab_red,
                                method="REML",
                                test="t")


ggplot(fauna_ab_red,aes(year.c,yi,size=1/vi))+
  geom_point(alpha=0.3)+
  theme_cowplot()
#some indication that effect sizes increase over time


#run an all-in publication bias test
abundance_red_all_in_bias_model<-rma.mv(yi,vi,
                                        mods = ~sei+
                                                year.c+
                                                perc_annual_dist-1, 
                                      random=~1|Site_ID/Study_ID,
                                      data=fauna_ab_red,
                                      method="REML",
                                      test="t")




##########################################
#analysis of changes in alpha diversity###
##########################################

fauna_div_m0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=fauna_div)#null model
fauna_div_m1 <-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, random=~1|Site_ID/Study_ID,data=fauna_div)#impact of drought vs increases

#calculate cook distances
cooks_div_0<-cooks.distance(fauna_div_m0)
cooks_div_1<-cooks.distance(fauna_div_m1)
#combine into one dataframe
c_dists_div<-data.frame(cooks_div_0,cooks_div_1)
#compare cook distances for different models
ggplot(c_dists_div,aes(cooks_div_0,cooks_div_1))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline()

#again there are some cooks distances for model 1 that are much larger than for the null model
#this means that we should filter out high cooks distances for model 1 and then rerun the model
fauna_div_filtered<- fauna_div %>%cbind(cooks_div_1) %>%filter(cooks_div_1 < 3.0*mean(cooks_div_1,na.rm=TRUE))
#this removed 2 comparisons

#rerun analysis of impact of decreases vs decreases
fauna_div_m1_filtered<-rma.mv(yi,vi,mods = ~factor(disturbance_type)-1, 
                             random=~1|Site_ID/Study_ID,data=fauna_div_filtered)
#this shows a 9% reduction with precipitation decreases and a 12% increase for precipitation increasess

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
fauna_ab_filtered<-fauna_ab_filtered%>%
  select(-cooks_ab_1)
fauna_div_filtered<-fauna_div_filtered%>%
  select(-cooks_div_1)
fauna_filtered<-rbind(fauna_ab_filtered,fauna_div_filtered)

#relabel disturbance types
fauna_filtered<-fauna_filtered%>%
  mutate(disturbance=if_else(disturbance_type=="drought","Precipitation\nreduction","Precipitation\nincrease"),
         outcome=if_else(broad_outcome=="abundance","Abundance","Alpha diversity"))

fauna_filtered$broad_outcome

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
                  data = fauna_ab_filtered, 
                  xlab = "Change in abundance relative\nto baseline (log response ratio)",
                  k = TRUE, g = TRUE,k.pos="left")+
  theme_cowplot()+
  scale_fill_manual(values = c("#fde725","#1f9e89"))+
  scale_color_manual(values = c("#fde725","#1f9e89"))+
  scale_x_discrete(labels=c("Precipitation\ndecrease", "Precipitation\n increase"))+
  annotate("text", x = 2.3, y = 2, label = "+49%")+
  annotate("text", x = 1.3, y = 2, label = "-49%")+
  theme(legend.position = "none",
        text = element_text(size=10),
        axis.text = element_text(size=10))

#diversity model
orchard_div_plot<-orchard_plot(fauna_div_m1_filtered,  group ='Site_ID',mod = 'disturbance_type',
                                data = fauna_div_filtered, 
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



