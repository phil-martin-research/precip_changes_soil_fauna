# ORCHARD PLOT FIGURES

rm(list = ls())

pacman::p_load(tidyverse,tidyverse,cowplot,metafor,orchaRd,ggbeeswarm,tidyr,insight,gt,gtExtras,webshot)


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
abundance<-read_csv("data/abundance_data.csv")
diversity<-read_csv("data/diversity_data.csv")

#---------------------------------------------------
#1. data analysis
#---------------------------------------------------

#################################
#analysis of change in abundance#
#################################

########################################
#following reductions in precipitation##
########################################

fauna_ab_red_m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_red)#null model
#35% reduction in abundance

#remove comparisons that fail geary's test
abundance_red_no_geary<-abundance_red%>%
  mutate(geary_test=if_else(is.na(geary_test),"Not needed",geary_test))%>%
  filter(geary_test!="fail")

fauna_ab_red_m0_no_geary<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_red_no_geary)#null model
#32% reduction in abundance

#calculate cook distances
cooks_ab_red_0<-cooks.distance(fauna_ab_red_m0)

#filter out highly influential comparisons
abundance_red_filtered<-abundance_red%>%cbind(cooks_ab_red_0)%>%filter(cooks_ab_red_0 < 3.0*mean(cooks_ab_red_0,na.rm=TRUE))
#this removed 9 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_ab_red_m0_filter<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_red_filtered)#null model
#this shows a 32% reduction with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
abundance_red_appraisal<-abundance_red%>%
                         filter(Validity!="Low validity")

fauna_ab_red_m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_red_appraisal)#model with no low validity studies
#studies with higher robustness tended to show greater reductions in abundance,
#however studies with higher validity had higher reductions in precipitation

#put all this information into a table
#info to include - estimate, se, p val, Q result, I squared

#create loop to do this
model_type<-c("Null model","Outliers removed","Low validity removed")
model_list<-list(fauna_ab_red_m0,fauna_ab_red_m0_filter,fauna_ab_red_m0_no_low)
ab_dec_sensitivity_summary<-data.frame()
for (i in 1:3){
  params<-broom::tidy(model_list[[i]])
  qe<-model_list[[i]]$QE
  qe_p<-model_list[[i]]$QEp
  I2<-I2_multi(model_list[[i]])
  sens_temp<-data.frame(disturbance="Precipitation reduction",outcome="Abundance",
                        model_type=model_type[i],params,qe,qe_p,I2)
  ab_dec_sensitivity_summary<-rbind(ab_dec_sensitivity_summary,sens_temp)
}

########################################
#following increases in precipitation##
########################################

fauna_ab_inc_m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_inc)#null model
#this shows no significant increase with an increase of 19%

#calculate cook distances
cooks_ab_inc_0<-cooks.distance(fauna_ab_inc_m0)

#filter out highly influential comparisons
abundance_inc_filtered<- abundance_inc %>%cbind(cooks_ab_inc_0) %>%filter(cooks_ab_inc_0 < 3.0*mean(cooks_ab_inc_0,na.rm=TRUE))
#this removed 7 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_ab_inc_m0_filter<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_inc_filtered)#null model
#this shows no significant increase with an increase of 12%

#run sensitivity analysis based on critical appraisal quality
abundance_inc_appraisal<-abundance_inc%>%
  filter(Validity!="Low validity")

fauna_ab_inc_m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_inc_appraisal)#model with no low validity studies
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
  sens_temp<-data.frame(disturbance="Precipitation increase",outcome="Abundance",
                        model_type=model_type[i],params,qe,qe_p,I2)
  sensitivity_summary_ab_inc<-rbind(sensitivity_summary_ab_inc,sens_temp)
}

##########################################
#analysis of changes in alpha diversity###
##########################################

########################################
#following reductions in precipitation##
########################################

fauna_div_red_m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_red)#null model

#calculate cook distances
cooks_div_red_0<-cooks.distance(fauna_div_red_m0)

#filter out highly influential comparisons
diversity_red_filtered<- diversity_red %>%cbind(cooks_div_red_0) %>%filter(cooks_div_red_0 < 3.0*mean(cooks_div_red_0,na.rm=TRUE))
#this removed 5 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_div_red_m0_filter<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_red_filtered)#null model
#this shows an 8% reduction of diversity with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
diversity_red_appraisal<-diversity_red%>%
  filter(Validity!="Low validity")

fauna_div_red_m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_red_appraisal)#model with no low validity studies
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
  sens_temp<-data.frame(disturbance="Precipitation reduction",outcome="Alpha diversity",
                        model_type=model_type[i],params,qe,qe_p,I2)
  div_red_sensitivity_summary<-rbind(div_red_sensitivity_summary,sens_temp)
}

########################################
#following increases in precipitation##
########################################

fauna_div_inc_m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_inc)#null model

#calculate cook distances
cooks_div_inc_0<-cooks.distance(fauna_div_inc_m0)

#filter out highly influential comparisons
diversity_inc_filtered<- diversity_inc %>%cbind(cooks_div_inc_0) %>%filter(cooks_div_inc_0 < 3.0*mean(cooks_div_inc_0,na.rm=TRUE))
#this removed 3 comparisons

#rerun analysis of impact of reductions in precipitation
fauna_div_inc_m0_filter<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_inc_filtered)#null model
#this shows a 4% reduction in diversity with precipitation decreases

#run sensitivity analysis based on critical appraisal quality
diversity_inc_appraisal<-diversity_inc%>%
  filter(Validity!="Low validity")

fauna_div_inc_m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_inc_appraisal)#model with no low validity studies
#similar reduction in effect size here

#studies with lower robustness tended to have a lower increase in precipitation
ggplot(diversity_inc,aes(x=Validity,y=perc_annual_dist))+
  geom_violin()

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
  sens_temp<-data.frame(disturbance="Precipitation increase",outcome="Alpha diversity",
                        model_type=model_type[i],params,qe,qe_p,I2)
  div_inc_sensitivity_summary<-rbind(div_inc_sensitivity_summary,sens_temp)
}


#combine tables of sensitivity analyses into one big table

sensitivity_results<-rbind(ab_dec_sensitivity_summary,sensitivity_summary_ab_inc,div_red_sensitivity_summary,div_inc_sensitivity_summary)


sensitivity_table <- sensitivity_results%>%
                     mutate(across(where(is.numeric), round, 3))%>%
                     gt()

sensitivity_table%>%gtsave("figures/for_paper/summary_sensitivity_table.docx")

#---------------------------------------------------
#3. Plot figures 
#---------------------------------------------------

#change in abundance

#make my own version of the plots
#bring together predictions from the different models
abun_inc_preds<-distinct(data.frame(predict(fauna_ab_red_m0)))
abun_dec_preds<-distinct(data.frame(predict(fauna_ab_inc_m0)))
div_inc_preds<-distinct(data.frame(predict(fauna_div_red_m0)))
div_dec_preds<-distinct(data.frame(predict(fauna_div_inc_m0)))
comb_preds<-rbind(abun_inc_preds,abun_dec_preds,div_inc_preds,div_dec_preds)
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

fauna_data<-rbind(abundance,diversity)

#relabel disturbance types
fauna_data<-fauna_data%>%
  mutate(disturbance=if_else(disturbance_type=="drought","Precipitation\nreduction","Precipitation\nincrease"),
         outcome=if_else(broad_outcome=="abundance","Abundance","Alpha diversity"))

#make my own orchard plot
ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=fauna_data,aes(x=lnrr_laj,y=disturbance,colour=outcome,group=outcome,size=1/v_lnrr_laj),
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
facet_summary_plot<-ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=subset(fauna_data,lnrr_laj<2.5),
                   aes(x=lnrr_laj,y=disturbance,colour=disturbance,group=outcome,size=1/v_lnrr_laj),
                   dodge.width = 0.01,alpha=0.5)+
  geom_errorbarh(data=comb_preds_2,aes(y=disturbance,xmin=pi.lb,xmax=pi.ub),
                 position=position_dodge(width=1),size=1,height=0,colour="black",alpha=0.8)+
  geom_errorbarh(data=comb_preds_2,aes(xmin=ci.lb,xmax=ci.ub,y=disturbance),
                 position=position_dodge(width=1),size=2,height=0,colour="black",alpha=0.8)+
  geom_point(data=comb_preds_2,aes(x=pred,y=disturbance,colour=disturbance,fill=disturbance),
             position=position_dodge(width=1),size=4,shape=21,colour="black")+
  theme_cowplot()+
  facet_wrap(~outcome,scales = "free_x")+
  scale_fill_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_color_manual("Outcome type",values = c("#fde725","#1f9e89"))+
  scale_size_continuous(range = c(1,10))+
  labs(y="Disturbance type",x="Soil and litter fauna relative to no disturbance (log response ratio)")+
  guides(size = "none")+
  theme(text=element_text(size=12),
        axis.text=element_text(size=10),
        legend.position = "bottom",
        legend.justification = "centre")

#create data for the labelling using number of studies and number of comparisons
studies_label<-fauna_data%>%
  group_by(disturbance,outcome)%>%
  summarise(k_count=length(lnrr_laj),
            studies=n_distinct(Study_ID),
            studies_annotation=paste("k = ",k_count," (",studies,")",sep=""))
studies_label$lnrr_laj<-c(-3,-0.6,-3,-0.6)

#create data for labelling using effect sizes in percentages etc
effect_size_label<-comb_preds_2%>%
            mutate(perc_change=round(perc_pred,0),
            change_label=if_else(perc_change<0,
                                 paste(perc_change,"%",sep=""),
                                 paste("+",perc_change,"%",sep="")))
effect_size_label$lnrr_laj<-c(1.3,1.3,0.4,0.4)

#add these to the plot
facet_summary_plot_with_label1<-facet_summary_plot+
  geom_text(data=studies_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=studies_annotation),
                hjust   = 0.3,
                vjust   = -4)


facet_summary_plot_with_label1+
  geom_text(data=effect_size_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=change_label),
            hjust   = 0,
            vjust   = -4)


ggsave("figures/for_paper/abun_div_summary_facet.pdf",width = 20,height = 14,units = "cm",dpi = 300)

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



