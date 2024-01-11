#script to draw conceptual figures

rm(list = ls())

#load packages
pacman::p_load(tidyverse,lemon,cowplot,patchwork,gghighlight)


#1 - differences between precipitation reduction and precipitation increase

#create data
precip_diff<-data.frame(dist_type=c("Precipitation\nreduction","Precipitation\nincrease"),
           es=c(-1,+1),
           lci=c(-1.5,0.5),
           uci=c(-0.5,1.5))

#plot this
precip_diff_plot<-ggplot(precip_diff,aes(x=es,y=dist_type,colour=dist_type))+
  geom_errorbarh(data=precip_diff,aes(xmin=lci,xmax=uci),colour="black",height=0)+
  geom_point(size=3,shape=19,aes(colour=dist_type))+
  geom_point(size=3,shape=21,colour="black")+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  labs(x="Effect size",y="")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fdf925"))
precip_diff_plot
ggsave("figures/for_paper/conceptual_diagram/disturbance_type.png",width = 5,height = 4,units = "cm",dpi=300)

#2 - change with precipitation

#create data
precip_change<-data.frame(precip_change=seq(-100,100,1))

precip_change_plot<-ggplot(precip_change,aes(precip_change,precip_change))+
  geom_line(linewidth=1.5,alpha=0.5)+
  theme_cowplot()+
  labs(y="Effect size",x="Change in precipitation")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)
precip_change_plot
ggsave("figures/for_paper/conceptual_diagram/precip_change.png",width = 5,height = 4,units = "cm",dpi=300)

#3 - modification of precipitation effect by microhabitat
precip_microhabitat_change<-data.frame(precip_change=rep(seq(-100,100,1),each=2),
                                       microhabitat=rep(c("Litter","Soil"),201))

#make litter more responsive to change than soils
precip_microhabitat_change$yi<-ifelse(precip_microhabitat_change$microhabitat=="Litter",
                                      precip_microhabitat_change$precip_change*1,
                                      precip_microhabitat_change$precip_change*0.5)

microhabitat_text<-data.frame(precip_change=c(-40,20,60),microhabitat=c("None","Litter","Soil"),
                              yi=c(100,100,100),label=c("Fauna habitat:","Litter","Soil"))

microhabitat_plot<-ggplot(precip_microhabitat_change,aes(precip_change,yi,colour=microhabitat))+
  geom_line(alpha=0.8,linewidth=1.5)+
  theme_cowplot()+
  labs(y="Effect size",x="Change in precipitation")+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  scale_x_continuous(limits = c(-100,105))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual(values=c("#70e48d","#d934a1"))+
  scale_fill_manual(values=c("#70e48d","white","#d934a1"))
microhabitat_plot
ggsave("figures/for_paper/conceptual_diagram/microhabitat_change.png",width = 5,height = 4,units = "cm",dpi=300)

#4 modification of precipitation effect by fauna size

#create data
precip_size_change<-data.frame(precip_change=rep(seq(-100,100,1),each=3),
                                       faunal_size=rep(c("Microfauna","Mesofauna","Macrofauna"),201))
precip_size_change<-precip_size_change%>%mutate(faunal_size=fct_relevel(faunal_size,"Microfauna","Mesofauna","Macrofauna"))

#make micro and mesofauna more responsive to change than macrofauna
precip_size_change$yi<-ifelse(precip_size_change$faunal_size=="Macrofauna",
                              precip_size_change$precip_change*0.1,
                              ifelse(precip_size_change$faunal_size=="Mesofauna",
                                     precip_size_change$precip_change*0.2,
                                     precip_size_change$precip_change*0.3))

#create text labels
size_text<-data.frame(precip_change=c(0,0,0),faunal_size=c("Microfauna","Mesofauna","Macrofauna"),
                              yi=c(50,50,50))
size_text<-size_text%>%mutate(faunal_size=fct_relevel(faunal_size,"Microfauna","Mesofauna","Macrofauna"))


microhabitat_text<-data.frame(precip_change=c(-40,20,60),microhabitat=c("None","Litter","Soil"),
                              yi=c(100,100,100),label=c("Fauna habitat:","Litter","Soil"))

#draw figure
precip_size_plot<-ggplot(precip_size_change,aes(precip_change,yi,colour=faunal_size))+
  geom_line(alpha=0.8,linewidth=1.5)+
  labs(y="Effect size",x="Change in precipitation")+
  scale_x_continuous(limits = c(-100,105))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  scale_color_manual(values = c("#3e728c","#df7f9a","#ffe64b"))
precip_size_plot

ggsave("figures/for_paper/conceptual_diagram/size_change.png",width = 10,height = 4,units = "cm",dpi=300)

#5 - exoskeleton alters impacts of precipitation changes
precip_exo_change<-data.frame(precip_change=rep(seq(-100,100,1),each=2),
                                       exoskeleton=rep(c("No exoskeleton","Exoskeleton"),201))

#make litter more responsive to change than soils
precip_exo_change$yi<-ifelse(precip_exo_change$exoskeleton=="No exoskeleton",
                             precip_exo_change$precip_change*1,
                             precip_exo_change$precip_change*0.2)

exo_text<-data.frame(precip_change=c(40,40),exoskeleton=c("No exoskeleton","Exoskeleton"),
                              yi=c(110,-20))

exo_plot<-ggplot(precip_exo_change,aes(precip_change,yi,colour=exoskeleton))+
  geom_line(alpha=0.8,linewidth=1.5)+
  geom_vline(xintercept=0,lty=2)+
  geom_hline(yintercept=0,lty=2)+
  theme_cowplot()+
  labs(y="Effect size",x="Change in precipitation")+
  scale_x_continuous(limits = c(-100,100))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual(values=c("#9a48ff","#ff9800"))
exo_plot
ggsave("figures/for_paper/conceptual_diagram/exoskeleton_change.png",width = 5,height = 4,units = "cm",dpi=300)

#6 put all plots together into one diagram

#do the same but using cowplot
cowplot_row_1<-plot_grid(NULL,precip_diff_plot,precip_change_plot,NULL,ncol=4,
                         labels=c("","(a)","(b)",""),label_size = 9,rel_widths = c(0.3,1.2,1,0.5))

cowplot_row_2<-plot_grid(microhabitat_plot,exo_plot,precip_size_plot,ncol=3,labels =c("(c)","(d)","(e)"),label_size=9)

cowplot_combined<-plot_grid(cowplot_row_1,cowplot_row_2,ncol=1)
cowplot_combined
save_plot("figures/for_paper/conceptual_diagram/conceptual_figure_multi.png",cowplot_combined,base_height = 10,base_width = 15,units="cm")  


#2 - conceptual figure for discussion

body_size<-data.frame(body_size=seq(0,100))

#add data on refuge ability
body_size$refuge<-1-body_size$body_size

#add data on physical resistance
body_size$resistance<-(body_size$body_size-50)^2

#add data on population change
body_size$pop_change<-(((body_size$body_size-50)^2)-2500)/2500

#plot graphs for drought

#refuge
refuge_plot<-ggplot(body_size,aes(body_size,refuge))+
  geom_line()+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="",y="Refuge ability")

#physical resistance
resistance_plot<-ggplot(body_size,aes(body_size,resistance))+
  geom_line()+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="",y="Physical resistance to drought")

#population change
population_plot_drought<-ggplot(body_size,aes(body_size,pop_change))+
  geom_line()+
  theme_cowplot()+
  geom_hline(yintercept = 0,lty=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="",y="Population change")+
  expand_limits(y=1)


#do the same but using cowplot
drought_row_1<-plot_grid(refuge_plot,resistance_plot,population_plot_drought,ncol=3,
                         labels=c("(a)","(b)","(c)"),label_size = 12)

ggsave("figures/for_paper/conceptual_diagram/drought_hypotheses.png",width = 12,height = 4,dpi = 300)
