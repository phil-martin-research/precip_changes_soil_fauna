#script to prepare data for different analyses:
#1- dataset for meta-analyses
#2 - spatial data for maps etc
#3 - dataset for rainfall bias etc

rm(list = ls())

library(tidyverse)
library(metafor)
library(tidyr)

#read in .csv file with soil fauna data
crit_appraisal<- read_csv("data/critical_appraisal_2023_05_29.csv")
sites<- read_csv("data/site_data_2023_05_29.csv")
fact_table<- read_csv("data/fact_table_2023_05_29.csv")
taxonomy<- read_csv("data/taxonomy.csv")
body_length<- read_csv("data/body_length.csv")
body_width<- read_csv("data/body_width.csv")


#---------------------------------------------------
#1 - format datasets for meta-analysis
#---------------------------------------------------

#clean dataset
fact_table <- fact_table %>%
  mutate(
    #remove disturbance strengths that represent >500% increases in precipitation
    perc_annual_dist = if_else(perc_annual_dist < 500, perc_annual_dist, NA_real_),
    perc_during_dist = if_else(perc_during_dist < 500, perc_during_dist, NA_real_),
    #change all variances that are 0 to NA
    control_var = if_else(control_var !=0, control_var, NA_real_),
    dist_var = if_else(dist_var !=0, dist_var, NA_real_),
    #alter format of precipitation change to represent continuous change from increase to decrease
    perc_annual_dist=if_else(perc_annual_dist>100, perc_annual_dist-100, perc_annual_dist*(-1)),
    perc_during_dist=if_else(perc_during_dist>100, perc_during_dist-100, perc_during_dist*(-1)),
    #convert SE to SD
    control_SD=if_else(var_type=="SE", control_var*sqrt(control_n), control_var),
    dist_SD=if_else(var_type=="SE",dist_var*sqrt(dist_n), dist_var)
  )


#join with data from sites and from critical appraisal
fact_table<-fact_table%>%
  left_join(sites,"Site_ID")%>%
  left_join(crit_appraisal,"Study_ID")%>%
  left_join(taxonomy,"Highest_taxonomic_resolution")%>%
  left_join(body_length, by=c('body_length'='Taxonomy'))%>%
  left_join(body_width, by=c('body_width'='Taxonomy'))


#remove columns that we don't use 
col_details<-data.frame(col_name=names(fact_table),
                        col_index=seq(1,116))

fact_table<-select(fact_table,-c(9,10,15,16,21,22,30,33:35,39:48,52:73,75:80,82:91,93:101,103:106,108:109,111:116))

#subset to remove data that represents a subset of other data
fact_table <- filter(fact_table, is_subset == FALSE)

#check to see if any of the means are equal to zero
#control group
fact_table%>%
  group_by(control_av)%>%
  summarise(length(control_av))
#disturbance group
fact_table%>%
  group_by(disturbance_av)%>%
  summarise(length(disturbance_av))

#there are only 19 data points where control or disturbance mean are equal to 0
#so we will exclude these
fact_table<-fact_table%>%
  filter(control_av>0&disturbance_av>0)

#extract year of study
fact_table$Study_ID<-ifelse(fact_table$Study_ID=="Aupic-Samain_2021a","Aupic_Samain_2021a",fact_table$Study_ID)
fact_table$study_year<-parse_number(fact_table$Study_ID,trim_ws = TRUE)

##############################################
#I need to sort this part out################
#############################################

#can I actually use this method to impute variance?

#impute missing SD values
alpha <- 0.05
z_crit <- qnorm(1 - alpha / 2)


fact_table<-fact_table%>%
  mutate(exact_p_val=as.numeric(exact_p_val),
         approx_p_value_numeric=ifelse(approx_p_val==">0.1",median(c(1,0.1)),NA),
         approx_p_value_numeric=ifelse(approx_p_val==">0.09",median(c(1,0.09)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.08",median(c(1,0.08)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.07",median(c(1,0.07)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.06",median(c(1,0.06)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.05",median(c(1,0.05)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.1",0.1,approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.5",0.05,approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.01",0.01,approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.001",0.001,approx_p_value_numeric))


#back calculate the pooled SD from means, group sizes, and the p-value based on this 
#code from wolfgang https://gist.github.com/wviechtb/3834b707caf50948f15c314d2a26a84c


fact_table<-fact_table%>%
  mutate(exact_tval=qt(exact_p_val/2, df=control_n+dist_n-2, lower.tail=FALSE),
         exact_pooled_SD=abs((disturbance_av - control_av) / (exact_tval * sqrt(1/dist_n + 1/control_n))),
         approx_tval=qt(approx_p_value_numeric/2, df=control_n+dist_n-2, lower.tail=FALSE),
         approx_pooled_SD=abs((disturbance_av - control_av) / (approx_tval * sqrt(1/dist_n + 1/control_n))),
         pooled_SD_to_use=if_else(!is.na(exact_pooled_SD),exact_pooled_SD,approx_pooled_SD),
         logRR_var=pooled_SD_to_use*((1/(dist_n*(disturbance_av^2)))+(1/(control_n*(control_av^2)))))



#use the recommended method of Nakagawa et al for calculation of variance



#impute values
#based on the median coefficient of variation
#and missing sample sizes based on median sample sizes
control_cv<-median(fact_table$control_SD/fact_table$control_av,na.rm = TRUE)
disturbance_cv<-median(fact_table$dist_SD/fact_table$disturbance_av,na.rm = TRUE)
med_control_n<-median(fact_table$control_n,na.rm = TRUE)
med_dist_n<-median(fact_table$dist_n,na.rm = TRUE)

fact_table<-fact_table%>%
  mutate(control_SD=if_else(is.na(control_SD),control_cv*control_av,control_SD),
         dist_SD=if_else(is.na(dist_SD),disturbance_cv*disturbance_av,dist_SD),
         control_n=if_else(is.na(control_n),med_control_n,control_n),
         dist_n=if_else(is.na(dist_n),med_dist_n,dist_n))

#calculate log response ratio
soil_fauna_rr<- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                       sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = fact_table)


ggplot(soil_fauna_rr,aes(vi,logRR_var))+
  geom_point()+
  geom_abline()+
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")

#add variable for the standard error of each effect size for publication bias analysis
soil_fauna_rr$sei <- sqrt(soil_fauna_rr$vi)

#calculate effective sample size for publication bias analysis
soil_fauna_rr$e_n<-with(soil_fauna_rr,(4*(control_n*dist_n)) / (control_n + dist_n))

#calculate the inverse of the "effective sample size" to account for unbalanced sampling for publication bias analysis
soil_fauna_rr$inv_n_tilda <-with(soil_fauna_rr, (control_n + dist_n)/(control_n*dist_n))
soil_fauna_rr$sqrt_inv_n_tilda <- with(soil_fauna_rr, sqrt(inv_n_tilda))


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

##################################################################
#2- format spatial data###########################################
##################################################################

#we want data with study details, locations, taxonomic groups/size classes, and disturbance types

