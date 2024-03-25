#this script produces an analysis of how different trophic groups of nematode are impacted by precipitation change

rm(list = ls())

pacman::p_load(tidyverse,metafor,cowplot,MuMIn,orchaRd)

#load data
all_fauna<-read.csv("data/fauna_data.csv")


##########################################
#tidy data################################
##########################################

nematode_trophic<-all_fauna%>%
  filter(Highest_taxonomic_resolution=="Nematoda",
         detailed_outcome=="abundance",
         use_for_trophic_analysis==TRUE,
         disturbance_type=="drought")


ggplot(nematode_trophic,aes(trophic_to_use,lnrr_laj))+
  geom_beeswarm()

# create a unit-level random effect to model residual variance in metafor
nematode_trophic$obsID <- 1:nrow(nematode_trophic)

#run model
nem_trophic_M1<-rma.mv(lnrr_laj,v_lnrr_laj,mods = ~trophic_to_use-1,random=~1|Study_ID/Site_ID/obsID,data=nematode_trophic)

#get model results
nem_results<-mod_results(nem_trophic_M1,mod="trophic_to_use",group="Study_ID")

#plot results
orchard_plot(nem_results,mod="trophic_to_use",group="Study_ID",xlab="Change in abundance (lnRR)")+
  theme(axis.text.y = element_text(angle = 0))
ggsave("figures/for_paper/nematode_for_reviewer.png",width = 18,height = 15,units="cm",dpi=300)
