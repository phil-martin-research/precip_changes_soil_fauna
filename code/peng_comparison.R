#figure for presentation

pacman::p_load(tidyverse,cowplot,metafor,orchaRd,ggbeeswarm,tidyr,insight,gt,gtExtras,webshot,scales,egg,lemon)


#make figure to compare our study results to those of Peng et al

peng_comp<-read.csv("data/peng_comparison.csv")

ggplot(peng_comp,aes(y=Disturbance,x=summary,xmin=lower,xmax=upper,colour=Study))+
  geom_vline(xintercept = 0,lty=2)+
  geom_point(position=position_dodge(width=0.5),size=3)+
  geom_errorbarh(position=position_dodge(width=0.5),height = .5,size=1)+
  facet_wrap(~Outcome,scales = "free_x")+
  xlab("Change relative to control (%)")+
  theme_cowplot()+
  scale_colour_viridis_d(begin = 0.5,end = 0.8)
ggsave("figures/for_presentations/peng_comparison.png",width = 25,height = 15,dpi = 300,units = "cm")  

#plot just the results of peng et al

peng_comp%>%
  filter(Study=="Peng et al")%>%
  ggplot(aes(y=Disturbance,x=summary,xmin=lower,xmax=upper))+
  geom_vline(xintercept = 0,lty=2)+
  geom_point(position=position_dodge(width=0.5),size=3)+
  geom_errorbarh(position=position_dodge(width=0.5),height = .5,size=1)+
  facet_rep_wrap(~Outcome,repeat.tick.labels = TRUE)+
  xlab("Change relative to control (%)")+
  theme_cowplot()+
  scale_colour_viridis_d(begin = 0.5,end = 0.8)

ggsave("figures/for_presentations/peng_results.png",width = 25,height = 10,dpi = 300,units = "cm")  
