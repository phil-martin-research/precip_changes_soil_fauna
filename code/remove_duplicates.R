rm(list = ls())


#code to find and remove duplicate articles from final search for work

library(synthesisr)
library(tidyverse)
library(janitor)

#locate bibliographic files
bibfiles<- list.files("data/searches/final_search/",
                       full.names = TRUE,recursive = FALSE,
                       pattern="*ris"
)


#read in refs
imported_files <- read_refs(
  filename = bibfiles,
  return_df = TRUE)

# first, we will remove articles that have identical titles
# this is a fairly conservative approach, so we will remove them without review
df <- deduplicate(
  imported_files,
  match_by = "title",
  method = "exact"
)

# then we will use string distance to identify likely duplicates
duplicates_string <- find_duplicates(
  df$title,
  method = "string_osa",
  to_lower = TRUE,
  rm_punctuation = TRUE,
  threshold = 7
)

#remove duplicates based on string distance
results <- extract_unique_references(df, duplicates_string)


#read in all references screened as part of systematic map
df_sm<- read_refs(
  filename = "data/searches/deduplicated_sm_searches.ris",
  return_df = TRUE)
#label studies as being from the systematic map
df_sm$search_source<-"systematic map"

#join systematic map references to those from systematic review and remove those already screened
sm_sr_join<-results%>%left_join(df_sm,"title")%>%
  filter(is.na(search_source))%>%
  select(1:23)%>%
  rename_with(~str_remove(., '.x'))


#write ris files
new_studies_ris<-write_refs(sm_sr_join, format="ris",file = "data/searches/final_search/clean_references/new_references.ris")
