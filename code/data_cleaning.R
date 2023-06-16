#script to prepare data for different analyses:
#1- dataset for meta-analyses
#2 - spatial data for maps etc
#3 - dataset for rainfall bias etc

rm(list = ls())

pacman::p_load(tidyverse,metafor,tidyr,here,patchwork,dplyr,raster)

#read in .csv file with soil fauna data
crit_appraisal<- read_csv("data/critical_appraisal_2023_05_29.csv")
sites<- read_csv("data/site_data_2023_06_07.csv")
fact_table<- read_csv("data/fact_table_2023_06_09.csv")
taxonomy<- read_csv("data/taxonomy.csv")
body_length<- read_csv("data/body_length.csv")
body_width<- read_csv("data/body_width.csv")


#---------------------------------------------------
#1 - format datasets for meta-analysis
#---------------------------------------------------

#join with data from sites and from critical appraisal
fact_table<-fact_table%>%
  left_join(sites,"Site_ID")%>%
  dplyr::select(-Study_ID.y)%>%
  rename(Study_ID=Study_ID.x)%>%
  left_join(crit_appraisal,by="Study_ID")

unique(fact_table$Highest_taxonomic_resolution)

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
    control_SD=ifelse(var_type=="SE", control_var*sqrt(control_n), control_var),
    dist_SD=ifelse(var_type=="SE",dist_var*sqrt(dist_n), dist_var),
    lat=if_else(is.na(Lat_dec_deg),Latitude_deg+(Latitude_min/60)+(Latitude_sec/3600),Lat_dec_deg),
    lon=if_else(is.na(Lon_dec_deg),Longitude_deg+(Longitude_min/60)+(Longitude_sec/3600),Lon_dec_deg),
    trophic_to_use=if_else(is.na(trophic_group),trophic_assigned,trophic_group),
    exoskeleton=if_else(Highest_taxonomic_resolution=="Nematoda"|
                        Highest_taxonomic_resolution=="Testate amoebae"|
                        Highest_taxonomic_resolution=="Tylenchidae"|
                        Highest_taxonomic_resolution=="Criconematidae"|
                        Highest_taxonomic_resolution=="Aphelenchoididae"|
                        Highest_taxonomic_resolution=="Aphelenchoididae"|
                        Highest_taxonomic_resolution=="Cephalobidae"|
                        Highest_taxonomic_resolution=="Plectidae"|
                        Highest_taxonomic_resolution=="Qudsianematidae",
                        "no","yes"))
#we have 661 rows here

#remove columns that we don't use 
col_details<-data.frame(col_name=names(fact_table),
                        col_index=seq(1,131))

fact_table<-dplyr::select(fact_table,-c(16,17,20,21,33,35:75,77:85,87:106,113:120))

#check to see if any of the means are equal to zero
#control group
fact_table%>%
  group_by(control_av)%>%
  summarise(length(control_av))
#disturbance group
fact_table%>%
  group_by(disturbance_av)%>%
  summarise(length(disturbance_av))

#there are only 10 data points where control or disturbance mean are equal to 0
#so we will exclude these - leaving us with 646 rows
fact_table<-fact_table%>%
  filter(control_av>0&disturbance_av>0)

#extract year of study
fact_table$study_year<-parse_number(fact_table$Study_ID,trim_ws = TRUE)

################################################
#2 - data imputation############################
################################################

#work out exact or approximate p values
fact_table<-fact_table%>%
  mutate(exact_p_val=as.numeric(exact_p_val),
         approx_p_value_numeric=ifelse(approx_p_val==">0.1",median(c(1,0.1)),NA),
         approx_p_value_numeric=ifelse(approx_p_val==">0.09",median(c(1,0.09)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.08",median(c(1,0.08)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.07",median(c(1,0.07)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.06",median(c(1,0.06)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val==">0.05",median(c(1,0.05)),approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.1",0.1,approx_p_value_numeric),
         approx_p_value_numeric=ifelse(approx_p_val=="<0.05",0.05,approx_p_value_numeric),
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
#based on scripts from https://alistairmcnairsenior.github.io/Miss_SD_Sim/

# Calculate CV on missing dataset. Note missing data will be ignored
fact_table <- fact_table %>%
  mutate(cv_Control = na_if(control_SD / control_av, Inf),
         cv_Experimental = na_if(dist_SD / disturbance_av, Inf))

# Function to calculate Geary's "number"
geary <- function(mean, sd, n){
  (1 / (sd / mean)) * ((4*n)^(3/2) / (1 + 4*n))
}

# Geary's test; assumption of normality assumed to be approximately correct when values are >= 3.
fact_table <- fact_table %>% 
  mutate(geary_control = geary(control_av, control_SD, control_n),
         geary_trt = geary(disturbance_av, dist_SD, dist_n),
         geary_test = ifelse(geary_control >= 3 & geary_trt >= 3, "pass", "fail"))
# How many fail?
geary_res <- fact_table %>% group_by(geary_test) %>% summarise(n = n()) %>%  data.frame()
#24 effect sizes fail representing around 5% of the data


# Calculate the average between study CV, which will replace missing values.
fact_table <- cv_avg(x = control_av, sd = control_SD,
                n = control_n, group = Study_ID, label = "1",
                data = fact_table)
fact_table <- cv_avg(x = disturbance_av, sd = dist_SD,
                n = dist_n, group = Study_ID,
                label = "2", data = fact_table)

# Use weighted mean CV in replacement for where CV's are missing. Otherwise, calculate CV^2 of data that is known.
fact_table <- fact_table %>%
  mutate(cv2_cont_new = if_else(is.na(cv_Control),      b_CV2_1, cv_Control^2),
         cv2_expt_new = if_else(is.na(cv_Experimental), b_CV2_2, cv_Experimental^2))

# Now calculate new yi and vi, called lnrr_laj & v_lnrr_laj, respectively. 
#This uses either the between individual CV^2 when missing or normal CV^2 when not missing.
fact_table <- fact_table %>%
  mutate(lnrr_laj = -lnrr_laj(m1 = control_av, m2 = disturbance_av,
                             cv1_2 = cv2_cont_new, cv2_2 = cv2_expt_new,
                             n1= control_n, n2 = dist_n),
         v_lnrr_laj = v_lnrr_laj(cv1_2 = cv2_cont_new, n1= control_n,
                                 cv2_2 = cv2_expt_new, n2 = dist_n))

# We need to exclude some missing data in the raw data set and data that is not defined on the ratio scale.
fact_table <- fact_table %>% filter(!is.infinite(lnrr_laj) & !is.na(lnrr_laj))

#calculate log response ratio
soil_fauna_rr<- escalc(m2i = control_av, m1i = disturbance_av, n2i = control_n, n1i = dist_n,
                       sd2i = control_SD, sd1i = dist_SD,  measure = "ROM", data = fact_table)


ggplot(soil_fauna_rr,aes(yi,lnrr_laj))+
  geom_point()+
  geom_abline()

soil_fauna_rr$rr_diff<-(abs(as.numeric(soil_fauna_rr$yi)-soil_fauna_rr$lnrr_laj))


ggplot(soil_fauna_rr,aes(rr_diff))+
  geom_histogram()+
  facet_wrap(~geary_test)

#the similarity of the two effect sizes is very high - need to work out where the deviations come from though
#it looks like they come from the comparisons that fail the geary test, we will test the impact of including these later


#add variable for the standard error of each effect size for publication bias analysis
soil_fauna_rr$sei <- sqrt(soil_fauna_rr$vi)

#calculate effective sample size for publication bias analysis
soil_fauna_rr$e_n<-with(soil_fauna_rr,(4*(control_n*dist_n)) / (control_n + dist_n))

#calculate the inverse of the "effective sample size" to account for unbalanced sampling for publication bias analysis
soil_fauna_rr$inv_n_tilda <-with(soil_fauna_rr, (control_n + dist_n)/(control_n*dist_n))
soil_fauna_rr$sqrt_inv_n_tilda <- with(soil_fauna_rr, sqrt(inv_n_tilda))


soil_fauna_rr<-as_tibble(soil_fauna_rr)

###################################################################
#3 - add aridity data##############################################
###################################################################


aridity_index<-raster("data/spatial_data/aridity/Global-AI_ET0_v3_annual/ai_v3_yr.tif")

#extract aridity data
soil_fauna_rr$aridity<-extract(aridity_index,cbind(soil_fauna_rr$lon,soil_fauna_rr$lat))
#convert to correct units
soil_fauna_rr$aridity<-soil_fauna_rr$aridity/10000

#remove columns that are not needed
col_details<-data.frame(col_name=names(soil_fauna_rr),
                        col_index=seq(1,75))

soil_fauna_rr<-dplyr::select(soil_fauna_rr,-c(7:11,21:25,38:41,43:44,50:60,62:65))


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

