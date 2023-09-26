#script to draw conceptual figures

rm(list = ls())

#load packages
pacman::p_load(tidyverse,lemon,cowplot,patchwork)


#1 - differences between precipitation reduction and precipitation increase

#create data
precip_diff<-data.frame(dist_type=c("Precipitation\nreduction","Precipitation\nincrease"),
           es=c(-1,+1),
           lci=c(-1.5,0.5),
           uci=c(-0.5,1.5))

#plot this
precip_diff_plot<-ggplot(precip_diff,aes(x=es,y=dist_type,colour=dist_type))+
  geom_pointrange(data=precip_diff,aes(xmin=lci,xmax=uci))+
  geom_point(size=3,shape=21,colour="black")+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()+
  labs(x="Effect size",y="")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual("Disturbance type",values = c("#1f9e89","#fde725"))
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
        axis.text = element_text(size=8))
precip_change_plot
ggsave("figures/for_paper/conceptual_diagram/precip_change.png",width = 5,height = 4,units = "cm",dpi=300)

#3 - modification of precipitation effect by microhabitat
precip_microhabitat_change<-data.frame(precip_change=rep(seq(-100,100,1),each=2),
                                       microhabitat=rep(c("Litter","Soil"),201))

#make litter more responsive to change than soils
precip_microhabitat_change$yi<-ifelse(precip_microhabitat_change$microhabitat=="Litter",
                                      precip_microhabitat_change$precip_change*1,
                                      precip_microhabitat_change$precip_change*0.5)

microhabitat_text<-data.frame(precip_change=c(90,90),microhabitat=c("Litter","Soil"),
                              yi=c(110,65),label=c("Litter","Soil"))

microhabitat_plot<-ggplot(precip_microhabitat_change,aes(precip_change,yi,colour=microhabitat))+
  geom_line(alpha=0.5,linewidth=1.5)+
  theme_cowplot()+
  labs(y="Effect size",x="Change in precipitation")+
  geom_text(data=microhabitat_text,aes(label=label),size=3)+
  scale_x_continuous(limits = c(-100,105))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual(values=c("#39b15f","#d934a1"))
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

#draw figure
precip_size_plot<-ggplot(precip_size_change,aes(precip_change,yi,colour=faunal_size))+
  geom_line(alpha=0.5,linewidth=1.5)+
  labs(y="Effect size",x="Change in precipitation")+
  geom_text(data=size_text,aes(label=faunal_size),size=3)+
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
  facet_rep_wrap(~faunal_size,repeat.tick.labels = TRUE)+
  scale_color_manual(values = c("#02475f","#c3386b","#e0b500"))
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
  geom_line(alpha=0.5,linewidth=1.5)+
  theme_cowplot()+
  labs(y="Effect size",x="Change in precipitation")+
  geom_text(data=exo_text,aes(label=exoskeleton),size=3)+
  scale_x_continuous(limits = c(-100,100))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        text=element_text(size=8),
        axis.text = element_text(size=8))+
  scale_color_manual(values=c("#39b15f","#d934a1"))
exo_plot
ggsave("figures/for_paper/conceptual_diagram/exoskeleton_change.png",width = 5,height = 4,units = "cm",dpi=300)

#6 put all plots together into one diagram
column_1<-plot_spacer()+precip_diff_plot+plot_spacer()+plot_layout(ncol=1)
column_2<-plot_spacer()+precip_change_plot+plot_spacer()+plot_layout(ncol=1)
column_3<-microhabitat_plot+precip_size_plot+exo_plot+plot_layout(ncol=1)

#horizontal
plot_spacer()+plot_spacer()+microhabitat_plot+
precip_diff_plot+precip_change_plot+precip_size_plot+
plot_spacer()+plot_spacer()+exo_plot

ggsave("figures/for_paper/conceptual_diagram/conceptual_figure_h.png",width = 20,height = 12,units="cm",dpi=300)


#vertical
plot_spacer()+precip_diff_plot+plot_spacer()+
plot_spacer()+precip_change_plot+plot_spacer()+
microhabitat_plot+precip_size_plot+exo_plot+
  plot_annotation(tag_levels = list(c('(a)', '(b)',"(c)","(d)","(e)")))+
  +plot_layout(widths = c(1, 1))


ggsave("figures/for_paper/conceptual_diagram/conceptual_figure_v.png",width = 15,height = 12,units="cm",dpi=300)

#or just as a multi-panel plot
precip_diff_plot+precip_change_plot+microhabitat_plot+exo_plot+precip_size_plot+
  plot_layout(ncol=2,widths=c(1,1,1,1,2))+
  plot_annotation(tag_levels = list(c('(a)', '(b)',"(c)","(d)","(e)")))
ggsave("figures/for_paper/conceptual_diagram/conceptual_figure_multi.png",width = 15,height = 15,units="cm",dpi=300)

#do the same but using cowplot
row_1_2<-plot_grid(precip_diff_plot,precip_change_plot,microhabitat_plot,exo_plot,
                   ncol=2,labels =c("(a)","(b)","(c)","(d)"),label_size=9)
multi_cowplot<-plot_grid(row_1_2,precip_size_plot,ncol=1,labels=c("","(e)"),rel_heights = c(2,1),label_size=9)

save_plot("figures/for_paper/conceptual_diagram/conceptual_figure_multi2.png",multi_cowplot,base_height = 15,base_width = 12,units="cm")  

plot_spacer()+precip_diff_plot+plot_spacer()+plot_layout(ncol=1)
column_2<-plot_spacer()+precip_change_plot+plot_spacer()+plot_layout(ncol=1)
column_3<-microhabitat_plot+precip_size_plot+exo_plot+plot_layout(ncol=1)
