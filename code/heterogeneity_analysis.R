# script to explore heterogeneity in the effects of precipitation changes on soil and litter fauna 

rm(list = ls())

pacman::p_load(tidyverse,metafor,cowplot,orchaRd,ggbeeswarm,tidyr,ggthemes,sp,broom,lemon)

#read in .csv files with soil fauna data
abundance<- read_csv("data/abundance_data.csv")
diversity<- read_csv("data/diversity_data.csv")


###############################################################################
#1 - data tidying##############################################################
###############################################################################

abundance%>%
  group_by(Validity,perc_annual_dist,Functional_group_size.y,aridity,exoskeleton)%>%
  summarise(size_count=length(above_below))%>%
  print(n=100)


#need to fix problem with validity

abundance$Validity<-ifelse(is.na(abundance$Validity),"Medium validity",abundance$Validity)

#complete cases for the variable about percentage change in precipitation and body size
abundance_complete<-abundance[complete.cases(abundance$perc_annual_dist,abundance$Functional_group_size.y,
                                             abundance$aridity,abundance$above_below,
                                             abundance$exoskeleton,abundance$Validity),]
#we lose 10 comparisons for abundance
diversity_complete<-diversity[complete.cases(diversity$perc_annual_dist,diversity$Functional_group_size),]
#we lose 7 comparisons for diversity


###############################################################################
#1 - models of heterogeneity in response of abundance of soil and litter fauna#
###############################################################################

#first a saturated model including all potential predictors
M_sat_abun<-rma.mv(lnrr_laj,v_lnrr_laj,mods =  ~aridity*perc_annual_dist+
                                                Functional_group_size.y*perc_annual_dist+
                                                above_below*perc_annual_dist+
                                                exoskeleton*perc_annual_dist+
                                                Validity*perc_annual_dist,
                   random=~1|Site_ID/Study_ID,data=abundance_complete)
#use cooks distance to identify influential points
cooks_sat_abun<-cooks.distance(M_sat_abun)
#filter out high cooks distances for saturated model and then run all models
abundance_filtered<- abundance_complete %>%cbind(cooks_sat_abun) %>%filter(cooks_sat_abun < 3.0*mean(cooks_sat_abun,na.rm=TRUE))
#this removes 11 comparisons - should we be doing this?

#test all models against each other
M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_complete)
M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist ,random=~1|Site_ID/Study_ID,data=abundance_complete)
M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*aridity,random=~1|Site_ID/Study_ID,data=abundance_complete)
M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size,random=~1|Site_ID/Study_ID,data=abundance_complete)
M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below,random=~1|Site_ID/Study_ID,data=abundance_complete)
M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size+above_below,random=~1|Site_ID/Study_ID,data=abundance_complete)


M4_no_coll<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*above_below,random=~1|Site_ID/Study_ID,data=abundance_complete_non_coll)

M4
M4_no_coll


vif(M1)
vif(M2)
vif(M3)
vif(M4)
vif(M5)

AIC(M1,M2,M3,M4,M5)

#centre perc_annual_dist and aridity
abundance_complete$perc_annual_dist_centred<-abundance_complete$perc_annual_dist-mean(abundance_complete$perc_annual_dist)
abundance_complete$aridity_centred<-abundance_complete$aridity-mean(abundance_complete$aridity,na.rm=TRUE)

#rerun models with new centred variables
M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=abundance_complete,method="ML")
M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred ,random=~1|Study_ID/Site_ID,data=abundance_complete,method="ML")
M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*aridity_centred,random=~1|Study_ID/Site_ID,data=abundance_complete,method="ML")
M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*Functional_group_size,random=~1|Study_ID/Site_ID,data=abundance_complete,method="ML")
M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*above_below,random=~1|Study_ID/Site_ID,data=abundance_complete,method="ML")
M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*Functional_group_size+I(perc_annual_dist_centred^2)*Functional_group_size,random=~1|Study_ID/Site_ID,data=abundance_complete,method="ML")



vif(M1)
vif(M2)
vif(M3)
vif(M4)
vif(M5)

ggplot(abundance_complete,aes(perc_annual_dist_centred,lnrr_laj,colour=above_below))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~above_below)


#calculate cook distances
cooks_M4<-cooks.distance(M4)

#filter out highly influential comparisons
abundance_complete_filtered<- abundance_complete %>%cbind(cooks_M4) %>%filter(cooks_M4 < 3.0*mean(cooks_M4,na.rm=TRUE))

M4_filtered<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*above_below,random=~1|Study_ID/Site_ID,data=abundance_complete_filtered,method="ML")


ggplot(abundance_complete_filtered,aes(perc_annual_dist_centred,lnrr_laj,colour=above_below))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~above_below)

AIC.rma(M0,M1,M2,M3,M4,M5)

#variance inflation factors all seem good now and M3 is the best model

#rerun the null and best models using REML
M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Study_ID/Site_ID,data=abundance_complete,method="REML")
M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*Functional_group_size-1,
           random=~1|Study_ID/Site_ID,data=abundance_complete,method="REML")


#calculate model fit
(M0$sigma2 - M3$sigma2) / M0$sigma2

#copy and save the model formula
M3_formula<-(~perc_annual_dist_centred*Functional_group_size-1)

#create dataframe with new data for predictions
new_data<-data.frame(expand.grid(
  perc_annual_dist_centred = seq(min(abundance_complete$perc_annual_dist_centred),max(abundance_complete$perc_annual_dist_centred),0.1),
  Functional_group_size=levels(as.factor(abundance_complete$Functional_group_size))))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M3_formula,data=new_data)

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(M3,newmods=predgrid))

#attach predictions to variables for plotting
new_data <- cbind(new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#plot data for model
#first subset to remove very large effect sizes
abundance_sub<-abundance_complete%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))

#remove predictions for microfauna and macrofauna where the changes in precipitation are outside the observed values
abundance_sub%>%
  group_by(Functional_group_size)%>%
  summarise(min_precip=min(perc_annual_dist_centred),
            max_precip=max(perc_annual_dist_centred))

new_data_micro<-subset(new_data,Functional_group_size=="microfauna"&perc_annual_dist_centred<=68.2)
new_data_meso<-subset(new_data,Functional_group_size!="microfauna")



#merge these together
new_data_merge<-rbind(new_data_micro,new_data_meso)
#create new variable for annual precipitation change
new_data_merge$perc_annual_dist<-new_data_merge$perc_annual_dist+mean(abundance_complete$perc_annual_dist)

#new version of size and annual change plot
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
ggplot(aes(x=perc_annual_dist,y=pred,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=ci.ub,ymin=ci.lb),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=pi.ub,ymin=pi.lb),colour=NA)+
  facet_rep_wrap(~Functional_group_size,repeat.tick.labels = TRUE)+
  geom_point(data=abundance_sub,aes(x=perc_annual_dist,y=lnrr_laj,size=1/v_lnrr_laj),alpha=0.25)+
  xlab("Change in annual precipitation (%)")+
  ylab("Change in soil fauna abundance (log ratio)")+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))+
  scale_fill_manual(values = c("#02475f","#c3386b","#e0b500"))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0,lty=2,alpha=0.3)+
  geom_vline(xintercept = 0,lty=2,alpha=0.3)
ggsave("figures/for_paper/abundance_precip_size.png",width = 20,height = 10,units = "cm",dpi = 300)

#same plot but with percentage change on the abundance axis
new_data_merge%>%
  mutate(Functional_group_size=fct_relevel(Functional_group_size,"microfauna","mesofauna","macrofauna"))%>%
  ggplot(aes(x=perc_annual_dist,y=(exp(pred)-1)*100,colour=Functional_group_size,fill=Functional_group_size))+
  geom_line()+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(ci.ub)-1)*100,ymin=(exp(ci.lb)-1)*100),colour=NA)+
  geom_ribbon(alpha=0.25,aes(ymax=(exp(pi.ub)-1)*100,ymin=(exp(pi.lb)-1)*100),colour=NA)+
  facet_wrap(~Functional_group_size,scales = "free")+
  geom_point(data=abundance_sub,aes(x=perc_annual_dist,y=(exp(lnrr_laj)-1)*100,size=1/v_lnrr_laj),alpha=0.25)+
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
abundance_red_all_in_bias_model<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~-1+sqrt_inv_n_tilda+
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

#first a saturated model including all potential predictors
M_sat_div<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~aridity+perc_annual_dist+I(perc_annual_dist^2)+Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_complete)
#use cooks distance to identify influential points
cooks_sat_div<-cooks.distance(M_sat_div)
#filter out high cooks distances for saturated model and then run all models
diversity_filtered<- diversity_complete %>%cbind(cooks_sat_div) %>%filter(cooks_sat_div < 3.0*mean(cooks_sat_div,na.rm=TRUE))
#this removes 3 comparisons

#test all models against each other
M0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M2<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M3<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~aridity*perc_annual_dist,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M4<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~aridity*Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M5<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist+Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M6<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M7<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M8<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M9<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist*aridity*Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")

#check to see which model is the best fit
AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9)
vif(M3)
vif(M6)
#the correlation is between the two is pretty weak though
ggplot(diversity_filtered,aes(perc_annual_dist,aridity))+
  geom_point()+
  geom_smooth(method = "lm")

#there are problems of multicolliniarity with models using the 'perc_annual_dist' and 'aridity' variables
#to attempt to deal with this, we will centre both variables my subtracting the mean from them
#magnitude of changes first
diversity_filtered$perc_annual_dist_centred<-diversity_filtered$perc_annual_dist-mean(diversity_filtered$perc_annual_dist)
#then aridity
diversity_filtered$aridity_centred<-diversity_filtered$aridity-mean(diversity_filtered$aridity)

ggplot(diversity_filtered,aes(perc_annual_dist_centred,aridity_centred))+
  geom_point()+
  geom_smooth(method = "lm")

#rerun all models and check variance inflation again
M0_centred<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M1_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M2_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M3_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~aridity_centred*perc_annual_dist_centred,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M4_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~aridity_centred*Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M5_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred+Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M6_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M7_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")
M8_centred<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~perc_annual_dist_centred*aridity_centred+Functional_group_size,random=~1|Site_ID/Study_ID,data=diversity_filtered,test="t",dfs="contain",method="ML")


AIC(M1_centred,M2_centred,M3_centred,M4_centred,M5_centred,M6_centred,M7_centred,M8_centred)

vif(M3_centred)
vif(M4_centred)
vif(M5_centred)
vif(M6_centred)
vif(M7_centred)
vif(M8_centred)

#model 3, an interaction between change in precipitation and aridity, seems to be the best

#copy and save the model formula
M3_formula<-(~1+aridity_centred*perc_annual_dist_centred)

#create dataframe with new data for predictions
summary(diversity_filtered$aridity_centred)
summary(diversity_filtered$perc_annual_dist_centred)

(1.13+0.7)/100
(85+90)/100

new_data<-data.frame(expand.grid(perc_annual_dist_centred = seq(-85,90,0.175),
                                 aridity_centred=seq(-1.13,0.7,0.00183)))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M3_formula,data=new_data)[,-1]

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(M3_centred,newmods=predgrid))

#attach predictions to variables for plotting
new_data <- cbind(new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#add original variables
new_data$perc_annual_dist<-new_data$perc_annual_dist_centred+mean(diversity_filtered$perc_annual_dist)
new_data$aridity<-new_data$aridity_centred+mean(diversity_filtered$aridity)

#create a convex hull of the change in precipitation and aridity data
ch<-chull(diversity_filtered$perc_annual_dist,diversity_filtered$aridity)
poly.df <- diversity_filtered[c(ch, ch[1]),]
poly <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(poly.df[,c(35,77)]))),1)))

# Create a SpatialPointsDataFrame with new.data:
new_data_poly <- new_data
coordinates(new_data_poly) <- ~perc_annual_dist+aridity

# Extract the points in new.data which are covered by the polygon:
new_data$inp <- over(new_data_poly, poly)
new_data <- new_data[complete.cases(new_data[,10]),]


#plot values as a heatmap
ggplot(new_data,aes(perc_annual_dist,aridity,fill=pred))+
  geom_raster()+
  scale_fill_gradient2_tableau(palette = "Orange-Blue-White Diverging",guide = "colorbar")+
  theme_cowplot()+
  scale_x_continuous(expand=expansion(0,0),limits = c(-100,85))+
  scale_y_continuous(expand=expansion(0,0))+
  labs(x="Change in annual precipitation (%)",y="Aridity index",fill="Change relative to\nbaseline (log response ratio)")+
  theme(text=element_text(size=12),
        axis.text=element_text(size=10),
        legend.title = element_text(size=10))+
  geom_jitter(data=diversity_filtered,aes(x=perc_annual_dist,y=aridity),inherit.aes = FALSE,shape=1,size=5)
ggsave("figures/for_paper/diversity_precip_arid_heatmap.png",width = 20,height = 14,units = "cm",dpi = 300)


#alternative version looking at effect sizes of different coefficients
tidy_M3<-tidy(M3_centred)

ggplot(tidy_M3,aes(x=estimate,xmax=estimate+(2*std.error),xmin=estimate-(2*std.error),y=term))+
  geom_point()+
  geom_errorbarh()+
  geom_vline(lty=2,xintercept = 0)+
  theme_cowplot()

#plot the interaction as different lines on a plot for different levels of aridity


#make predictions for aridity at 0.5, 1.0, 1.5, and 2.0
aridity_centred_groups<-c(0.5,1.0,1.5,2.0)-mean(diversity_filtered$aridity)

diversity_filtered%>%
  group_by(aridity_centred,perc_annual_dist_centred)%>%
  summarise(data_points=length(lnrr_laj))%>%
  ggplot(aes(perc_annual_dist_centred,aridity_centred,colour=data_points))+
  geom_point()


#create new dataset for predictions
new_data<-data.frame(perc_annual_dist_centred=c(seq(-80,50,1),seq(-20,75,1)),aridity_centred=c(rep(-0.75,times=131),rep(0.5,96)))

#create a model matrix and remove the intercept
predgrid<-model.matrix(M3_formula,data=new_data)[,-1]

#predict onto the new model matrix
mypreds<-data.frame(predict.rma(M3_centred,newmods=predgrid))

#attach predictions to variables for plotting
new_data <- cbind(new_data, mypreds[c("pred", "ci.lb", "ci.ub", "pi.lb", "pi.ub")])

#add original variables
new_data$perc_annual_dist<-new_data$perc_annual_dist_centred+mean(diversity_filtered$perc_annual_dist)
new_data$aridity<-new_data$aridity_centred+mean(diversity_filtered$aridity)


#plot the results of the predictions
ggplot(new_data,aes(perc_annual_dist,y=pred,ymin=ci.lb,ymax=ci.ub,colour=as.factor(aridity),fill=as.factor(aridity)))+
  geom_line(size=1)+
  geom_ribbon(alpha=0.5,colour=NA)+
  theme_cowplot()+
  labs(x="Change in annual precipitation (%)",
       y="Change in soil fauna alpha diversity\n(log response ratio)",
       fill="Aridity",colour="Aridity")+
  scale_color_manual(values = c("#1f78b4","#b2df8a"))+
  scale_fill_manual(values = c("#1f78b4","#b2df8a"))+
  facet_wrap(~aridity)+
  theme(legend.position = "none")
ggsave("figures/for_paper/diversity_precip_arid_lines.png",width = 20,height = 12,units = "cm",dpi = 300)
#not very convinced by this result and figure

#could divide between humid and not humid

diversity_filtered%>%
  mutate(arid_group=ifelse(aridity>0.65,"Humid","Not humid"))%>%
  ggplot(aes(x=perc_annual_dist))+
  geom_histogram()+
  facet_wrap(~arid_group)

diversity_filtered%>%
  mutate(arid_group=ifelse(aridity>0.65,"Humid","Not humid"))%>%
  ggplot(aes(x=perc_annual_dist,y=lnrr_laj))+
  geom_point(shape=1)+
  facet_wrap(~arid_group)+
  geom_smooth()
  


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


