#this script produces a series of descriptive statistics


pacman::p_load(tidyverse,cowplot,metafor)

#load data

abundance<-read.csv("data/abundance_data.csv")
fauna_data<-read.csv("data/fauna_data.csv")
sites<-read.csv("data/sites.csv")

#####################################################
#1 - how much data has information on size of fauna##
#####################################################

abundance_complete_width<-abundance[complete.cases(abundance$mean_width,abundance$perc_annual_dist,abundance$Functional_group_size.y,abundance$above_below,
                                                   abundance$exoskeleton,abundance$obsID,abundance$year.c),]

#plot histogram of this data
ggplot(abundance_complete_width,aes(x=mean_width))+
  geom_histogram()+
  scale_x_continuous(trans="log",breaks = c(0.1,0.5,1,2,5,10,20))+
  scale_y_continuous(expand = c(0, 0))+
  theme_cowplot()+
  labs(x="Mean body width (mm)",y="No. of effect sizes")
ggsave("figures/for_paper/body_width_histogram.png",width = 15,height = 10,units = "cm",dpi = 300)

######################################################
#2 - what season is sampling done in?#################
######################################################

unique(outcomes$sampling_season)

fauna_data%>%
  filter(detailed_outcome=="abundance"|detailed_outcome=="taxonomic richness"|detailed_outcome=="shannon wiener")%>%
  filter(!is.na(sampling_season))%>%
  filter(sampling_season!="")%>%
  group_by(sampling_season)%>%
  summarise(no_es=length(disturbance_av))%>%
  mutate(sampling_season=fct_relevel(sampling_season,"Spring","Summer","Autumn","Winter"))%>%
ggplot(aes(sampling_season,no_es,fill=sampling_season))+
  geom_bar(stat = "identity")+
  theme_cowplot()+
  scale_fill_manual(values = c("#c7f7a4",
                               "#fdee55",
                               "#f39e3e",
                               "#70afdf"))+
  theme(legend.position = "none")+
  labs(x="Sampling season",y="No. of effect sizes")
ggsave("figures/for_paper/sampling_season.png",width = 15,height = 10,units = "cm",dpi = 300)

#########################################################
#3 - Study length########################################
#########################################################

ggplot(fauna_data,aes(time_after_dist_start/365))+
  geom_histogram()+
  theme_cowplot()+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  labs(x="Time after precipitation change",y="No. of effect sizes")

total_es<-length(fauna_data$lnrr_laj)
fauna_data%>%
  filter(!is.na(time_after_dist_start))%>%
  mutate(time_after_dist=round(time_after_dist_start/365,0))%>%
  group_by(time_after_dist)%>%
  summarise(no_es=length(lnrr_laj),
            perc_es=(length(lnrr_laj)/total_es)*100)

############################################################
#4- country#################################################
############################################################

sites%>%
  group_by(Country)%>%
  summarise(site_count=length(Site_ID))

sites%>%
  summarise(site_count=length(Site_ID))


############################################################
#5-taxa sampled#############################################
############################################################

total_es<-length(fauna_data$lnrr_laj)

fauna_summmary<-fauna_data%>%
  group_by(Highest_taxonomic_resolution)%>%
  summarise(no_studies=length(lnrr_laj),
            perc_studies=(length(lnrr_laj)/total_es)*100)%>%
  arrange(desc(no_studies))%>%
  print(n=40)

fauna_summmary%>%
  mutate(other_fauna=if_else(Highest_taxonomic_resolution=="Collembola"|
                             Highest_taxonomic_resolution=="Oribatida"|
                             Highest_taxonomic_resolution=="Mesostigmata"|
                             Highest_taxonomic_resolution=="Acari"|
                             Highest_taxonomic_resolution=="Nematoda"|
                             Highest_taxonomic_resolution=="Diplopoda"|
                             Highest_taxonomic_resolution=="Prostigmata",
                             "No","Yes"))%>%
  group_by(other_fauna)%>%
  summarise(no_studies=sum(no_studies),
            perc_studies=sum(perc_studies),
            no_taxa=length(Highest_taxonomic_resolution))
  
(60+35+34+11)/429


############################################################
#6 - metric type############################################
############################################################

total_es<-length(fauna_data$lnrr_laj)

fauna_data%>%
  group_by(detailed_outcome)%>%
  summarise(no_es=length(lnrr_laj),
            perc_es=length(lnrr_laj)/total_es)

###############################################################
#7 -study type#################################################
###############################################################

study_length<-length(sites$exp_obs)

sites%>%
  group_by(exp_obs)%>%
  summarise(no_es=length(exp_obs),
            perc_es=length(exp_obs)/study_length)
  
  

