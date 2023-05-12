#In this script we test for publication bias following the recommendations of Nakagawa et al (2021) https://doi.org/10.1111/2041-210X.13724

rm(list = ls())

library(tidyverse)
library(cowplot)
library(metafor)
library(orchaRd)
library(ggbeeswarm)
library(tidyr)

#read in .csv files with soil fauna data
abundance<- read_csv("data/abundance_data.csv")
diversity<- read_csv("data/diversity_data.csv")
abundance_red<- read_csv("data/abundance_red_data.csv")
abundance_inc<- read_csv("data/abundance_inc_data.csv")
diversity_red<- read_csv("data/diversity_red_data.csv")
diversity_inc<- read_csv("data/diversity_inc_data.csv")

#############################################################################################
#1 - tests of publication bias for abundance under precipitation reduction###################
#############################################################################################


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
       ylab = "Effective sample size",
       xlab = "Effect size (lnRR)") 
#this does't look too bad


#use multilevel models to examine publication biases

# create a unit-level random effect to model residual variance in metafor
fauna_ab_red$obsID <- 1:nrow(fauna_ab_red)

#run an intercept-only model of impacts of precipitation reduction
fauna_ab_m0_reduction<-rma.mv(yi,vi,
                              random=list(~1|Site_ID/Study_ID,
                                          ~1|obsID),
                              method="REML",
                              test="t",
                              data=fauna_ab_red)

###############################################
#meta-regression with effective sampling size##
###############################################

#calculating the inverse of the "effective sample size" to account for unbalanced sampling

fauna_ab_red$inv_n_tilda <-  with(fauna_ab_red, (control_n + dist_n)/(control_n*dist_n))
fauna_ab_red$sqrt_inv_n_tilda <-  with(fauna_ab_red, sqrt(inv_n_tilda))

# Application of Equation 27 Nakagawa et al

fauna_ab_reduction_srin<-rma.mv(yi,vi,
                                mods=~1+sqrt_inv_n_tilda,
                                random=list(~1|Site_ID/Study_ID,
                                            ~1|obsID),
                                method="REML",
                                test="t",
                                data=fauna_ab_red)
#there is little evidence of a small-study effect, i.e. that effect sizes with larger uncertainty tend to be larger

###############################################
#meta-regression with year of publication######
###############################################

# mean-centering year of publication to help with interpretation
fauna_ab_red$year.c <- as.vector(scale(fauna_ab_red$study_year, scale = F))

#meta-regression of effect of publication year
abundance_red_year_bias_model<-rma.mv(yi,vi,mods = ~1+year.c, 
                                      random=list(~1|Site_ID/Study_ID,
                                                  ~1|obsID),
                                      data=fauna_ab_red,
                                      method="REML",
                                      test="t")


ggplot(fauna_ab_red,aes(year.c,yi,size=1/vi))+
  geom_point(alpha=0.3)+
  theme_cowplot()
#some indication that effect sizes increase over time


#run an all-in publication bias test
abundance_red_all_in_bias_model<-rma.mv(yi,vi,
                                        mods = ~-1+
                                          sqrt_inv_n_tilda+
                                          year.c+
                                          perc_annual_dist, 
                                        random=list(~1|Site_ID/Study_ID,
                                                    ~1|obsID),
                                        data=fauna_ab_red,
                                        method="REML",
                                        test="t")
#a slight suggestion a small-study effect

#plot this result
preds_abundance_red_all_in_bias_model<-data.frame(predict(abundance_red_all_in_bias_model,
                                                          newmods=cbind(seq(min(fauna_ab_red$sqrt_inv_n_tilda),
                                                                            max(fauna_ab_red$sqrt_inv_n_tilda),
                                                                            length.out=142),
                                                                        c(0),c(0)),addx=TRUE))

ggplot(preds_abundance_red_all_in_bias_model,aes(X.sqrt_inv_n_tilda,pred,ymin=ci.lb,ymax=ci.ub))+
  geom_line()+
  geom_ribbon(alpha=0.2)+
  geom_point(data=fauna_ab_red,aes(x=sqrt_inv_n_tilda,y=yi),alpha=0.2,inherit.aes = FALSE)+
  theme_cowplot()+
  labs(x="square root of inverse of effective sample size",y="effect size(lnRR)")
#impact is relatively small

