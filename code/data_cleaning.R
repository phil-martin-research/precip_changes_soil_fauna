#script to prepare data for different analyses:
#1- dataset for meta-analyses
#2 - spatial data for maps etc
#3 - dataset for rainfall bias etc

#first run file weighted_CV_functions.R

pacman::p_load(tidyverse,metafor,tidyr,here,patchwork,dplyr,raster,ggthemes,lemon)

#read in .csv file with soil fauna data
crit_appraisal<- read_csv("data/critical_appraisal.csv")
sites<- read_csv("data/sites.csv")
fact_table<- read_csv("data/outcomes.csv")
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
  left_join(crit_appraisal,by="Study_ID")%>%
  left_join(taxonomy,by="Highest_taxonomic_resolution")

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
                        "no","yes"),
    trap=if_else(sampling_method=="pitfall traps","Trap","Other"))
#we have 751 rows here

#remove columns that we don't use 
col_details<-data.frame(col_name=names(fact_table),
                        col_index=seq(1,89))

fact_table<-dplyr::select(fact_table,-c(18,19,22,23,33,36:45,47:68,70:78))

#check to see if any of the means are equal to zero
#control group
fact_table%>%
  group_by(control_av)%>%
  summarise(length(control_av))
#disturbance group
fact_table%>%
  group_by(disturbance_av)%>%
  summarise(length(disturbance_av))

#there are only 29 data points where control or disturbance mean are equal to 0
#so we will exclude these - leaving us with 725 rows
fact_table<-fact_table%>%
  filter(control_av>0&disturbance_av>0)

#extract year of study
fact_table$study_year<-parse_number(fact_table$Study_ID,trim_ws = TRUE)

#check the number of rows for which the sample size is missing
#and for which measure of variability is missing
fact_table%>%
  group_by(use_for_first_analysis)%>%
  summarise(perc_n_missing=(sum(is.na(control_n)&is.na(dist_n))/length(disturbance_av))*100,
            perc_var_missing=(sum(is.na(control_var)&is.na(dist_var))/length(disturbance_av))*100)

################################################
#2 - data imputation############################
################################################

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
#104 effect sizes fail representing around 15% of the data


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

#the similarity of the two effect sizes is very high
#it looks like the differences come from the comparisons that fail the geary test, we will test the impact of including these later


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
soil_fauna_rr$aridity<-raster::extract(aridity_index,cbind(soil_fauna_rr$lon,soil_fauna_rr$lat))
#convert to correct units
soil_fauna_rr$aridity<-soil_fauna_rr$aridity/10000

#remove columns that are not needed
col_details<-data.frame(col_name=names(soil_fauna_rr),
                        col_index=seq(1,63))

soil_fauna_rr<-dplyr::select(soil_fauna_rr,-c(9:12,23:27,50:53))

#check to see which outcomes are most commonly reported
soil_fauna_rr%>%
  group_by(detailed_outcome,disturbance_type)%>%
  summarise(no_es=length(lnrr_laj),no_studies=length(unique(Study_ID)))
#it only makes sense to analyse data on abundance, taxonomic richness, and shannon wiener

#subset dataset to get variables of interest for first set of analyses
#all relevant data
fauna_all <- filter(soil_fauna_rr, detailed_outcome=="abundance"|detailed_outcome=="taxonomic richness"|detailed_outcome=="shannon wiener")

# All abundance 
fauna_ab <- filter(soil_fauna_rr, broad_outcome == 'abundance',use_for_first_analysis==TRUE)
# Subset to give only studies of each precipitation change and outcome combination
fauna_ab_red <- filter(fauna_ab, disturbance_type == 'drought',use_for_first_analysis==TRUE)
fauna_ab_inc <- filter(fauna_ab, disturbance_type == 'precip_inc',use_for_first_analysis==TRUE)


# All taxonomic richness
fauna_richness <- filter(soil_fauna_rr, detailed_outcome == 'taxonomic richness',use_for_first_analysis==TRUE)
# Subset to give only studies of each precipitation change and outcome combination
fauna_rich_red <- filter(fauna_richness, disturbance_type == 'drought',use_for_first_analysis==TRUE)
fauna_rich_inc <- filter(fauna_richness, disturbance_type == 'precip_inc',use_for_first_analysis==TRUE)

# All Shannon wiener
fauna_shannon <- filter(soil_fauna_rr, detailed_outcome == 'shannon wiener',use_for_first_analysis==TRUE)
# Subset to give only studies of each precipitation change and outcome combination
fauna_shannon_red <- filter(fauna_shannon, disturbance_type == 'drought',use_for_first_analysis==TRUE)
fauna_shannon_inc <- filter(fauna_shannon, disturbance_type == 'precip_inc',use_for_first_analysis==TRUE)

#subset data to get abundance data to use for trophic analysis
# All abundance 
fauna_ab_trophic <- filter(soil_fauna_rr, broad_outcome == 'abundance',use_for_trophic_analysis==TRUE)


# Save all these files
write.csv(fauna_all, "data/fauna_data.csv")
write.csv(fauna_ab, "data/abundance_data.csv")
write.csv(fauna_richness, "data/richness_data.csv")
write.csv(fauna_shannon, "data/shannon_data.csv")
write.csv(fauna_ab_red, "data/abundance_red_data.csv")
write.csv(fauna_ab_inc, "data/abundance_inc_data.csv")
write.csv(fauna_rich_red, "data/richness_red_data.csv")
write.csv(fauna_rich_inc, "data/richness_inc_data.csv")
write.csv(fauna_shannon_red, "data/shannon_red_data.csv")
write.csv(fauna_shannon_inc, "data/shannon_inc_data.csv")
write.csv(fauna_ab_trophic, "data/trophic_abundance_data.csv")

##################################################################
#2- format spatial data###########################################
##################################################################

#we want data with study details, locations, taxonomic groups/size classes, and disturbance types
soil_fauna_spatial_data<-soil_fauna_rr%>%filter(use_for_first_analysis==TRUE)%>%#keep only data we use for analyses
                          dplyr::select(Study_ID,Site_ID,disturbance_type,Functional_group_size.y,perc_annual_dist,Country,lat,lon) #select only columns we want
                          
#save this as a .csv
write.csv(soil_fauna_spatial_data,"data/fauna_spatial_data.csv")



