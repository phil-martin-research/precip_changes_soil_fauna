#this script is produces summary meta-analyses and associated sensitivity tests
#for changes in abundance, taxonomic richness, and shannon diversity of soil and litter fauna
#as a result of increases and decreases in precipitation

rm(list = ls())


#load packages
pacman::p_load(tidyverse,tidyverse,cowplot,metafor,orchaRd,ggbeeswarm,tidyr,insight,gt,gtExtras,webshot,scales,egg,lemon)

######################################
#function to calculate I2 for multilevel models
I2_multi<-function(model){
  W <- diag(1/model$vi)
  X <- model.matrix(model)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
}

#####################################

#read in .csv files with soil fauna data
abundance<-read.csv("data/abundance_data.csv")
richness<-read.csv("data/richness_data.csv")
shannon<-read.csv("data/shannon_data.csv")
abundance_red<- read.csv("data/abundance_red_data.csv")
abundance_inc<- read.csv("data/abundance_inc_data.csv")
richness_red<- read.csv("data/richness_red_data.csv")
richness_inc<- read.csv("data/richness_inc_data.csv")
shannnon_red<- read.csv("data/shannon_red_data.csv")
shannnon_inc<- read.csv("data/shannon_inc_data.csv")

######################################################
#1. meta-analyses#####################################
######################################################

fauna_list<-list(abundance_red,abundance_inc,richness_red,richness_inc,shannnon_red,shannnon_inc)
outcomes<-c("Abundance","Abundance","Taxonomic richness","Taxonomic richness","Shannon-Wiener index","Shannon-Wiener index")
disturbances <- rep(c("Precipitation\nreduction", "Precipitation\nincrease"), times = 3)
sensitivity_summary<-data.frame()
prediction_summary<-data.frame()
#loop through all the different stages for the meta-analyses
for (i in 1:length(fauna_list)){
  #subset list to give relevant dataframe
  temp_df<-fauna_list[[i]]
  #run null model
  m0<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=temp_df)
  #remove comparisons that fail geary's test
  no_geary<-temp_df%>%
    mutate(geary_test=if_else(is.na(geary_test),"Not needed",geary_test))%>%
    filter(geary_test!="fail")
  #null model excluding the studies that fail Geary's test
  m0_no_geary<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=no_geary)
  #remove studies that have low validity as assessed by critical appraisal
  temp_df_appraisal<-temp_df%>%
    filter(Validity!="Low validity")
  #run model with no low validity studies
  m0_no_low<-rma.mv(lnrr_laj,v_lnrr_laj,random=~1|Site_ID/Study_ID,data=temp_df_appraisal)
  #put all this information into a table
  #info to include - estimate, se, p val, Q result, I squared
  #create loop to do this
  model_type<-c("Null model","Failed Geary test","Low validity removed")
  model_list<-list(m0,m0_no_geary,m0_no_low)
  for (y in 1:3){
    params<-broom::tidy(model_list[[y]])
    qe<-model_list[[y]]$QE
    qe_p<-model_list[[y]]$QEp
    I2<-I2_multi(model_list[[y]])
    k<-model_list[[y]]$k.all
    sens_temp<-data.frame(disturbance=disturbances[i],outcome=outcomes[i],
                          model_type=model_type[y],params,k,qe,qe_p,I2)
    sensitivity_summary<-rbind(sensitivity_summary,sens_temp)
  }
  #create predictions based on models
  temp_pred<-distinct(data.frame(predict(m0)))
  temp_pred2<-data.frame(disturbance=disturbances[i],detailed_outcome=outcomes[i],temp_pred)
  prediction_summary<-rbind(prediction_summary,temp_pred2)
}
  
#combine tables of sensitivity analyses into one big table
sensitivity_table <- sensitivity_summary%>%
                     mutate(across(where(is.numeric), round, 3))%>%
                     mutate(I2=round(I2,0))%>%
                     gt()

#export this to a word file
sensitivity_table%>%gtsave("figures/for_paper/summary_sensitivity_table.docx")


#convert predictions to percentages
prediction_summary_perc<-prediction_summary%>%
  mutate(perc_pred=(exp(pred)-1)*100,
         per_ci.lb=(exp(ci.lb)-1)*100,
         per_ci.ub=(exp(ci.ub)-1)*100,
         per_pi.lb=(exp(pi.lb)-1)*100,
         per_pi.ub=(exp(pi.ub)-1)*100)

#change names of different disturbances and outcomes

#---------------------------------------------------
#3. Plot figures 
#---------------------------------------------------

#produce my own version of an orchard plot

#first produce a plot for abundance

#Organise data into one dataset

fauna_data<-rbind(abundance,richness,shannon)

#relabel disturbance types
fauna_data<-fauna_data%>%
  mutate(disturbance=if_else(disturbance_type=="drought","Precipitation\nreduction","Precipitation\nincrease"),
         detailed_outcome=if_else(detailed_outcome=="abundance","Abundance",
                         if_else(detailed_outcome=="taxonomic richness","Taxonomic richness","Shannon-Wiener index")))%>%
  filter(lnrr_laj>-3)

#facetted version of the figure
  facet_summary_plot<-ggplot()+
  geom_vline(xintercept = 0,lty=2,size=1)+
  geom_quasirandom(data=subset(fauna_data,lnrr_laj<2.5),
                   aes(x=lnrr_laj,y=disturbance,colour=disturbance,group=detailed_outcome,size=1/v_lnrr_laj),
                   dodge.width = 0.01,alpha=0.5)+
  geom_errorbarh(data=prediction_summary_perc,aes(y=disturbance,xmin=pi.lb,xmax=pi.ub),
                 position=position_dodge(width=1),linewidth=1,height=0,colour="black",alpha=0.8)+
  geom_errorbarh(data=prediction_summary_perc,aes(xmin=ci.lb,xmax=ci.ub,y=disturbance),
                 position=position_dodge(width=1),linewidth=2,height=0,colour="black",alpha=0.8)+
  geom_point(data=prediction_summary_perc,aes(x=pred,y=disturbance,colour=disturbance,fill=disturbance),
             position=position_dodge(width=1),size=4,shape=21,colour="black")+
  theme_cowplot()+
  facet_rep_wrap(~detailed_outcome,scales = "free_x",repeat.tick.labels = TRUE,ncol=1)+
  scale_fill_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))+
  scale_size_continuous(range = c(1,10),transform = "sqrt")+
  labs(y="Disturbance type",x="Change in soil and litter fauna\noutcome (log response ratio)")+
  guides(size = "none")+
  theme(text=element_text(size=12),
        axis.text=element_text(size=10),
        legend.position = "none",
        legend.justification = "centre",
        strip.background =element_rect(fill="white",color = "black", linewidth = 1),
        strip.text = element_text(face = "bold",margin = unit(rep(2, 4), "pt")),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#create data for the labelling using number of studies and number of comparisons
studies_label<-fauna_data%>%
  group_by(disturbance,detailed_outcome)%>%
  summarise(k_count=length(lnrr_laj),
            studies=n_distinct(Study_ID),
            studies_annotation=paste("k = ",k_count," (",studies,")",sep=""))
studies_label$lnrr_laj<-c(-1.5,-1,-1,-1.5,-1,-1)

#create data for labelling using effect sizes in percentages etc
effect_size_label<-prediction_summary_perc%>%
            mutate(perc_change=round(perc_pred,0),
            change_label=if_else(perc_change<0,
                                 paste(perc_change,"%",sep=""),
                                 paste("+",perc_change,"%",sep="")),
            change_label=if_else(sign(per_ci.lb)==sign(per_ci.ub),paste(change_label,"*",sep = ""),change_label))
effect_size_label$lnrr_laj<-c(1.3,1.3,0.3,0.3,0.05,0.05)

#add these to the plot
facet_summary_plot_with_label1<-facet_summary_plot+
  geom_text(data=studies_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=studies_annotation),
                hjust   = 0.5,
                vjust   = -2.2)


facet_summary_plot_with_label2<-facet_summary_plot_with_label1+
  geom_text(data=effect_size_label,
            aes(x=lnrr_laj,
                y=disturbance,
                label=change_label),
            hjust   = 0,
            vjust   = -2.2)

#add labels for figures
facet_summary_plot_with_label3<-tag_facet(facet_summary_plot_with_label2)+
  theme(strip.text = element_text(face = "bold",margin = unit(rep(4, 4), "pt")),
        strip.background = element_rect(fill="white",color = "black", size = 1))
facet_summary_plot_with_label3

ggsave("figures/for_paper/abun_div_summary_facet.png",facet_summary_plot_with_label3,width = 16,height = 17,units = "cm",dpi = 300)
