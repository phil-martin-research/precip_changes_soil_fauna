#code for drawing map showing location of study sites

#load packages
pacman::p_load(tidyverse,ggthemes,lemon,raster,sf,sp,rgdal,dismo,statip,cowplot,plotrix,plotbiomes,ggnewscale,ggmagnify)

#load data
spatial_data<-read.csv("data/fauna_spatial_data.csv")

#####################################
#1 - map and biome plot##############
#####################################

####################################
#1.1 - map##########################
####################################


#summarise data to give numbers of comparisons from each site

spatial_data_unique<-spatial_data%>%mutate(precip_dec=if_else(disturbance_type=="drought",TRUE,FALSE),
                                           precip_inc=if_else(disturbance_type=="precip_inc",TRUE,FALSE))%>%
                                           group_by(lat,lon)%>%
                                           summarise(dec_count=sum(precip_dec),inc_count=sum(precip_inc))%>%
                                           mutate(dist_types=if_else(dec_count&inc_count>0,"Both increases\nand decreases",if_else(dec_count>0,
                                                 "Precipitation\ndecreases","Precipitation\nincreases")))%>%
                                           mutate(dist_types=fct_relevel(dist_types),"Precipitation\ndecreases")%>%
                                           group_by(lat,lon,dist_types)%>%
                                           summarise(total_comp=dec_count+inc_count)

#different size for different numbers of comparisons
#group sites into unique sites

#make global map
world_map <- map_data("world")
site_map<-ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray",colour="dark grey",linewidth=0.05)+
  theme_bw()+
  geom_jitter(data=spatial_data_unique,aes(x=lon,y=lat,shape=dist_types,size=total_comp,fill=dist_types),
             group=NA,alpha=0.7,stroke=0.2,width = 1)+
  coord_equal(ylim = c(-55,80))+
  theme(axis.title = element_text(colour="white"),
        axis.text = element_text(colour="white"),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.margin=margin(t=-30),
        text = element_text(size = 9),
        plot.margin = margin(0, 0, 0, -1,unit="cm"))+
  scale_y_continuous(expand = expansion(0,5))+
  scale_x_continuous(expand = expansion(0,0))+
  scale_shape_manual("Disturbance type",values = c(24,25,23))+
  scale_size_continuous(range = c(1,4),guide = 'none')+
  scale_fill_manual("Disturbance type",values = c("#7142ff","#ff412c","#e32eee"),guide="legend")+
  guides(fill = guide_legend(override.aes = list(size = 4)))

site_map

###########################################
#1.2 - biome plot##########################
###########################################


#now create plot for biomes
#turn coordinates into spatial points dataframe
xy <-dplyr::select(spatial_data_unique,lon,lat)


spdf <- SpatialPointsDataFrame(coords = xy, data = spatial_data_unique,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# ===== Prepare raster stack
# Read temperature and precipitation as raster stack.
# Low resolution raster datasets come with 'plotbiomes' package.
path <- system.file("extdata", "temp_pp.tif", package = "plotbiomes")
temp_pp <- raster::stack(path)
names(temp_pp) <- c("temperature", "precipitation")

# ===== Extract from raster stack
# Extract temperature and precipitation values from the raster datasets
extractions <- raster::extract(temp_pp, spdf, df = TRUE)
# Adjust temperature values to normal scale because WorldClim temperature data
# has a scale factor of 10 (integer storage for saving space).
extractions$temperature <- extractions$temperature/10

#join the dataframe of site characteristics to the extracted data
unique_site_climate<-cbind(spatial_data_unique,extractions)


#plot biome plot
biome_plot <- ggplot() +
  # add biome polygons
  geom_polygon(data = Whittaker_biomes,
               aes(x    = temp_c,
                   y    = precp_cm*10,
                   fill = biome),
               # adjust polygon borders
               colour = "gray98",
               linewidth   = 1)+
  # fill the polygons with predefined colors
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)+
  new_scale("fill")+
  #add site data
  geom_point(data = unique_site_climate, 
             aes(x = temperature, 
                 y = precipitation,
                 shape=dist_types,
                 size=total_comp,
                 fill=dist_types),
                 alpha  = 0.8)+
  ylab("Mean annual precipitation (mm)")+
  xlab("Mean annual temperature (°C)")+
  theme_bw()+
  theme(text = element_text(size = 10),
        axis.text = element_text(size=10),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=8),
        plot.margin = margin(0.3, 0, 0, 0,unit="cm"))+
  scale_y_continuous(expand = expansion(0,0))+
  scale_x_continuous(expand = expansion(0,0))+
  scale_shape_manual("Disturbance types",values = c(24,25,23),guide="none")+
  scale_size_continuous(range = c(1,4),guide = 'none')+
  scale_fill_manual("Disturbance types",values = c("#7142ff","#ff412c","#e32eee"),guide="none")
biome_plot

#combine the two plots
biome_combined<-plot_grid(biome_plot,NULL,
                          rel_widths = c(0.8,0.2))

#combined plots
map_biome<-plot_grid(site_map,NULL,biome_combined,
          labels=c("(a)","","(b)"),
          rel_heights = c(1.2,-0.2,1),
          nrow=3)

#save the figure
save_plot("figures/Figure_5.pdf",map_biome,base_height = 15,base_width = 18,dpi=300,units="cm")
