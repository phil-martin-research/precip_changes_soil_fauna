# Plot the relationship between abundance and strength of disturbance for both 
# drought and precipitation increase (all) including new body size data

#---------------------------------------------------
rm(list = ls())

library(ggpubr)
library(cowplot)
library(readxl)
library(ggplot2)
library(metafor)
library(dplyr)
library(tidyverse)
library(readr)
library(forcats)
library(raster)
library(maptools)

install.packages("vctrs")

Size_df <- read_csv("~/Google Drive/My Drive/holisoils_precip_inc_dec/data/joined_nosubset_env_length.csv")

#---------------------------------------------------
#adapt dataset for analysis
#---------------------------------------------------

#remove anomalies - outlying disturbance strengths 
Size_df$perc_annual_dist <- ifelse(Size_df$perc_annual_dist<500, Size_df$perc_annual_dist,  NA)
Size_df$perc_during_dist <- ifelse(Size_df$perc_during_dist<500, Size_df$perc_during_dist,  NA)
Size_df$av_length <- ifelse(Size_df$av_length<300, Size_df$av_length,  NA)

#change variances from 0 to NA
Size_df$control_var <- ifelse(Size_df$control_var !=0, Size_df$control_var,  NA)
Size_df$dist_var <- ifelse(Size_df$dist_var !=0, Size_df$dist_var, NA)

# change percentage to be continuous from decrease to increase 
Size_df$perc_annual_dist <- ifelse(Size_df$perc_annual_dist>100, Size_df$perc_annual_dist-100, Size_df$perc_annual_dist*(-1))
Size_df$perc_during_dist <- ifelse(Size_df$perc_during_dist>100, Size_df$perc_during_dist-100, Size_df$perc_during_dist*(-1))

#convert SE to SD
Size_df$control_SD<-ifelse(Size_df$var_type=="SE", Size_df$control_var*sqrt(Size_df$control_n), Size_df$control_var)
Size_df$dist_SD<-ifelse(Size_df$var_type=="SE",Size_df$dist_var*sqrt(Size_df$dist_n), Size_df$dist_var)


#---------------------------------------------------
#1. data exploration and effect sizes
#---------------------------------------------------
#log response ratio
Size_df_rr <- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                     sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = Size_df)

#impute missing variances with median variance 
Size_df_rr$vi<-ifelse(is.na(Size_df_rr$vi),median(Size_df_rr$vi, na.rm = TRUE), Size_df_rr $vi)

#subset Size_df_rr dataset to set variables of interest
size_ab_rr <- Size_df_rr %>%filter(broad_outcome == 'abundance')


#fit all models to same data: complete cases of variables and yi and filter cooks
#cuts the obs from 251 to 180 
complete_size_ab_rr <- size_ab_rr[complete.cases(size_ab_rr[ , c("Functional_group_size","av_length","av_width","perc_during_dist","perc_annual_dist","time_after_dist_start", "yi")]),]

complete_M0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=complete_size_ab_rr)
cooks<- cooks.distance(complete_M0)
size_ab <- complete_size_ab_rr %>%cbind(cooks) %>%filter(cooks < 3.0*mean(cooks))

# models 

names(size_ab)

ggplot(size_ab,aes(av_length,yi))+
  geom_point()
ggplot(size_ab,aes(av_width,yi))+
  geom_point()


M0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=size_ab)
M1<-rma.mv(yi,vi,mods = ~perc_annual_dist ,random=~1|Site_ID/Study_ID,data=size_ab)
M2<-rma.mv(yi,vi,mods = ~Functional_group_size ,random=~1|Site_ID/Study_ID,data=size_ab)
M3<-rma.mv(yi,vi,mods = ~perc_annual_dist+Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)
M4<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=size_ab)
M5<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=size_ab)
M6<-rma.mv(yi,vi,mods = ~perc_annual_dist*aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)
M7<-rma.mv(yi,vi,mods = ~av_width,random=~1|Site_ID/Study_ID,data=size_ab)
M8<-rma.mv(yi,vi,mods = ~perc_annual_dist*av_width+I(perc_annual_dist^2)*av_width-1,random=~1|Site_ID/Study_ID,data=size_ab)
M9<-rma.mv(yi,vi,mods = ~perc_annual_dist*av_width-1,random=~1|Site_ID/Study_ID,data=size_ab)


AIC_size_ab <- AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,
                       correct = TRUE)
AIC_size_ab<-AIC_size_ab[order(AIC_size_ab$AICc),]


#models inc. aridity 

M0<-rma.mv(yi,vi,random=~1|Site_ID/Study_ID,data=size_ab)
M1<-rma.mv(yi,vi,mods = ~perc_annual_dist ,random=~1|Site_ID/Study_ID,data=size_ab)
M2<-rma.mv(yi,vi,mods = ~Functional_group_size ,random=~1|Site_ID/Study_ID,data=size_ab)
M3<-rma.mv(yi,vi,mods = ~perc_annual_dist+Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)
M4<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=size_ab)
M5<-rma.mv(yi,vi,mods = ~perc_annual_dist*Functional_group_size+I(perc_annual_dist^2)*Functional_group_size-1,random=~1|Site_ID/Study_ID,data=size_ab)
M6<-rma.mv(yi,vi,mods = ~perc_annual_dist*aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)
M7<-rma.mv(yi,vi,mods = ~av_width,random=~1|Site_ID/Study_ID,data=size_ab)
M8<-rma.mv(yi,vi,mods = ~perc_annual_dist*av_width+I(perc_annual_dist^2)*av_width-1,random=~1|Site_ID/Study_ID,data=size_ab)
M9<-rma.mv(yi,vi,mods = ~perc_annual_dist*av_width-1,random=~1|Site_ID/Study_ID,data=size_ab)
M10<-rma.mv(yi,vi,mods = ~aridity+perc_annual_dist,random=~1|Site_ID/Study_ID,data=size_ab)
M11<-rma.mv(yi,vi,mods = ~aridity*perc_annual_dist,random=~1|Site_ID/Study_ID,data=size_ab)
M12<-rma.mv(yi,vi,mods = ~aridity+Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)
M13<-rma.mv(yi,vi,mods = ~aridity*Functional_group_size,random=~1|Site_ID/Study_ID,data=size_ab)

AIC_size_ab <- AIC.rma(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,
                       correct = TRUE)
AIC_size_ab<-AIC_size_ab[order(AIC_size_ab$AICc),]