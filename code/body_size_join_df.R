
#===============================================
#joining datasets together to include the body width and body size instead of highest taxonomic res
#===============================================
rm(list = ls())

#load libraries  ###########
library(tidyverse)
library(xlsx)

#---------------------------------------------------
# join data frame
#---------------------------------------------------

#load df

#load df

Env_df <- readxl::read_excel("~/Google Drive/My Drive/holisoils_precip_inc_dec/data/ALL_env_df_body_length.xlsx", sheet = 1)
body_length <- readxl::read_excel("~/Google Drive/My Drive/holisoils_precip_inc_dec/data/ALL_env_df_body_length.xlsx", sheet = 2)
body_width <- readxl::read_excel("~/Google Drive/My Drive/holisoils_precip_inc_dec/data/ALL_env_df_body_length.xlsx", sheet = 3)
taxa <- readxl::read_excel("~/Google Drive/My Drive/holisoils_precip_inc_dec/data/ALL_env_df_body_length.xlsx", sheet = 4)

#join

fact_taxa_df<-left_join(Env_df,taxa,by="Highest_taxonomic_resolution")
length_df<- left_join(fact_taxa_df, body_length, by=c('body_length'='Taxonomy'))
width_df<- left_join(length_df, body_width, by=c('body_width'='Taxonomy'))
join_all_na<- data.frame(lapply(width_df, function(x) {gsub("NA", "", x)}))

#write csv files
write.csv(join_all_na, 'data/join_Env_length_df.xlsx')

# subset == False data 
join_nosubset <- join_all_na %>%
  filter(is_subset =='FALSE')
write.csv(join_nosubset, "~/Google Drive/My Drive/holisoils_precip_inc_dec/data/study_data_for_analysis_26_04.csv") 


