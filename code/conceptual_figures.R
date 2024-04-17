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

#5 - exoskeleton alters impacts of precipitation changes
precip_exo_change<-data.frame(precip_change=rep(seq(-100,100,1),each=2),
                                       exoskeleton=rep(c("No exoskeleton","Exoskeleton"),201))

#make organisms with no exoskeleton more responsive to change than those with exoskeleton
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

#6 - aridity alters impacts of precipitation changes
precip_arid_change<-data.frame(precip_change=rep(seq(-100,100,1),each=2),
                               aridity=rep(c("Moist","Arid"),201))

#make organisms with no exoskeleton more responsive to change than those with exoskeleton
precip_arid_change$yi<-ifelse(precip_arid_change$aridity=="Moist",
                              precip_arid_change$precip_change*1.1+(-0.007*precip_arid_change$precip_change^2),
                              precip_arid_change$precip_change*0.6+(0.002*precip_arid_change$precip_change^2))

arid_text<-data.frame(precip_change=c(40,40),aridity=c("Arid","Moist"),
                     yi=c(110,-20))

arid_plot<-ggplot(precip_arid_change,aes(precip_change,yi,colour=aridity))+
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
  scale_color_manual(values=c("#E09C00","#3C19CA"))
arid_plot

#6 put all plots together into one diagram

#do the same but using cowplot
cowplot_row_1<-plot_grid(precip_diff_plot,precip_change_plot,microhabitat_plot,ncol=3,
                         labels=c("(a)","(b)","(c)"),label_size = 9)

cowplot_row_2<-plot_grid(exo_plot,precip_size_plot,arid_plot,ncol=3,labels =c("(d)","(e)","(f)"),label_size=9)

cowplot_combined<-plot_grid(cowplot_row_1,cowplot_row_2,ncol=1)
cowplot_combined
save_plot("figures/figure_1.png",cowplot_combined,base_height = 10,base_width = 15,units="cm")  
