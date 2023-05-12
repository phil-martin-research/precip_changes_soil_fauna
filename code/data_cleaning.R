#script to prepare data for different analyses

rm(list = ls())

library(tidyverse)
library(metafor)
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

#add variable for the standard error of each effect size for publication bias analysis
soil_fauna_rr$sei <- sqrt(soil_fauna_rr$vi)

#calculate effective sample size for publication bias analysis
soil_fauna_rr$e_n<-with(soil_fauna_rr,(4*(control_n*dist_n)) / (control_n + dist_n))

#calculate the inverse of the "effective sample size" to account for unbalanced sampling for publication bias analysis
soil_fauna_rr$inv_n_tilda <-  with(soil_fauna_rr, (control_n + dist_n)/(control_n*dist_n))
soil_fauna_rr$sqrt_inv_n_tilda <-  with(soil_fauna_rr, sqrt(inv_n_tilda))

#subset dataset to get variables of interest
# All abundance 
fauna_ab <- filter(soil_fauna_rr, broad_outcome == 'abundance')
# All diversity 
fauna_div <- filter(soil_fauna_rr, broad_outcome == 'alpha diversity')
# Subset to give only studies of each precipitation change and outcome combination
fauna_ab_red <- filter(fauna_ab, disturbance_type == 'drought')
fauna_ab_inc <- filter(fauna_ab, disturbance_type == 'precip_inc')
fauna_div_red <- filter(fauna_div, disturbance_type == 'drought')
fauna_div_inc <- filter(fauna_div, disturbance_type == 'precip_inc')

# Save all these files
write.csv(fauna_ab, "data/abundance_data.csv")
write.csv(fauna_div, "data/diversity_data.csv")
write.csv(fauna_ab_red, "data/abundance_red_data.csv")
write.csv(fauna_ab_inc, "data/abundance_inc_data.csv")
write.csv(fauna_div_red, "data/diversity_red_data.csv")
write.csv(fauna_div_inc, "data/diversity_inc_data.csv")
