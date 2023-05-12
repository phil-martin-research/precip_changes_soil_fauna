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
funnel(abundance_red$yi,
       abundance_red$vi,
       yaxis="seinv",
       ylab="Precision (1/SE)",
       xlab = "Effect size (lnRR)")
#this generally looks fine, although there is some bias towards negative values

#plot funnel plot with effective sample size
funnel(abundance_red$yi, abundance_red$vi, ni = abundance_red$e_n, yaxis="ni",
       ylab = "Effective sample size",
       xlab = "Effect size (lnRR)") 
#this doesn't look too bad

#multilevel models to examine publication biases

# create a unit-level random effect to model residual variance in metafor
abundance_red$obsID <- 1:nrow(abundance_red)

#meta-regression with effective sampling size using Equation 27 Nakagawa et al

fauna_ab_reduction_srin<-rma.mv(yi,vi,
                                mods=~1+sqrt_inv_n_tilda,
                                random=list(~1|Site_ID/Study_ID,
                                            ~1|obsID),
                                method="REML",
                                test="t",
                                data=abundance_red)
#there is little evidence of a small-study effect, i.e. that effect sizes with larger uncertainty tend to be larger

#meta-regression with year of publication

# mean-centering year of publication to help with interpretation
abundance_red$year.c <- as.vector(scale(abundance_red$study_year, scale = F))

#meta-regression of effect of publication year
abundance_red_year_bias_model<-rma.mv(yi,vi,mods = ~1+year.c, 
                                      random=list(~1|Site_ID/Study_ID,
                                                  ~1|obsID),
                                      data=abundance_red,
                                      method="REML",
                                      test="t")
#some indication that effect sizes increase over time


#run an all-in publication bias test
abundance_red_all_in_bias_model<-rma.mv(yi,vi,
                                        mods = ~-1+
                                          sqrt_inv_n_tilda+
                                          year.c+
                                          perc_annual_dist, 
                                          random=list(~1|Site_ID/Study_ID,
                                                    ~1|obsID),
                                          data=abundance_red,
                                          method="REML",
                                          test="t")
#a slight suggestion a small-study effect

#plot this result
preds_abundance_red_all_in_bias_model<-data.frame(predict(abundance_red_all_in_bias_model,
                                                          newmods=cbind(seq(min(abundance_red$sqrt_inv_n_tilda),
                                                                            max(abundance_red$sqrt_inv_n_tilda),
                                                                            length.out=142),
                                                                        c(0),c(0)),addx=TRUE))

ggplot(preds_abundance_red_all_in_bias_model,aes(X.sqrt_inv_n_tilda,pred,ymin=ci.lb,ymax=ci.ub))+
  geom_line()+
  geom_ribbon(alpha=0.2)+
  geom_point(data=abundance_red,aes(x=sqrt_inv_n_tilda,y=yi),alpha=0.2,inherit.aes = FALSE)+
  theme_cowplot()+
  labs(x="square root of inverse of effective sample size",y="effect size(lnRR)")
#impact is relatively small


#need to repeat this for:
# 1 - abundance and precipitation increase
# 2 - diversity and precipitation decreases
# 3 - diversity and precipitation increases
# 4 - abundance with increases and decreases combined
# 5 - diversity with increases and decreases combined


#############################################################################
#2 - tests of publication biases for all abundance data######################
#############################################################################

# create a unit-level random effect to model residual variance in metafor
abundance$obsID <- 1:nrow(abundance)

# mean-centering year of publication to help with interpretation
abundance$year.c <- as.vector(scale(abundance$study_year, scale = F))

#run an all-in publication bias test
abundance_all_in_bias_model<-rma.mv(yi,vi,
                                        mods = ~-1+
                                          sqrt_inv_n_tilda+
                                          year.c+
                                          perc_annual_dist, 
                                        random=list(~1|Site_ID/Study_ID,
                                                    ~1|obsID),
                                        data=abundance,
                                        method="REML",
                                        test="t")

abundance_model<-rma.mv(yi,vi,mods = ~perc_annual_dist, 
                                    random=list(~1|Site_ID/Study_ID,
                                                ~1|obsID),
                                    data=abundance,
                                    method="REML",
                                    test="t")

#evidence of impact of small-study effect - i.e. that effect sizes with larger uncertainty tend to be larger
#however, the impact of this on the slope related to change in precipitation is very small


#############################################################################
#3 - tests of publication biases for all diversity data######################
#############################################################################

# create a unit-level random effect to model residual variance in metafor
diversity$obsID <- 1:nrow(diversity)

# mean-centering year of publication to help with interpretation
diversity$year.c <- as.vector(scale(diversity$study_year, scale = F))

#run an all-in publication bias test
diversity_all_in_bias_model<-rma.mv(yi,vi,
                                    mods = ~-1+
                                      sqrt_inv_n_tilda+
                                      year.c+
                                      perc_annual_dist, 
                                    random=list(~1|Site_ID/Study_ID,
                                                ~1|obsID),
                                    data=diversity,
                                    method="REML",
                                    test="t")
#little evidence of impact of small-study effect or of a decline effect
