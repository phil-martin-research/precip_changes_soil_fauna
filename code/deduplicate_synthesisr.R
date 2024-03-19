#deduplicate and clear references for systematic map
library(synthesisr)
library(stringr)
library(tidyverse)

#tidy Google scholar references and export as RIS file
gs_scraped<-read.csv("data/full_search/scholar_scraping/scholar_scraped.csv")

gs_cleaned<-data.frame(
            label=NA,
            date_generated=as.Date("05/11/21"),
            source_type=NA,
            author=gsub(".*>","",gs_scraped$AU)[1],
            year=gs_scraped$PY,
            title=gsub(".*>","",gs_scraped$TI),
            journal=NA,
            volume=NA,
            issue=NA,
            start_page=NA,
            doi=gs_scraped$DO,
            issn=NA,
            url=gs_scraped$UR,
            id=paste(gsub(".*? ", "", str_match(gs_scraped$AU, ">\\s*(.*?)\\s*;")[,2]),gs_scraped$PY,sep = ""))

ris_out <- write_refs(gs_cleaned, format = "ris", file = "data/full_search/gs_scraped.ris")


#locate bibliographic files
bibfiles <- list.files("data/full_search_2/",
  full.names = TRUE,recursive = FALSE,
  pattern="*ris"
)

print(bibfiles)

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

#create small example to test the thresholds for removal of duplicates

samp_df<-sample_n(df,2000)

matches_theshold<-data.frame()
manual_checks<-data.frame()
#loop to test the effect of increasing the threshold
for (n in 1:20){
  duplicates_string_n <- find_duplicates(
    samp_df$title,
    method = "string_osa",
    to_lower = TRUE,
    rm_punctuation = TRUE,
    threshold = n
  )
  manual_checks_n<- review_duplicates(samp_df$title, duplicates_string_n)
  manual_checks_n$n<-n
  manual_checks<-rbind(manual_checks,manual_checks_n)
  matches_theshold_n<-data.frame(theshold=n,
             n_matches=length(unique(manual_checks_n$matches)))
  matches_theshold<-rbind(matches_theshold,matches_theshold_n)
}


ggplot(matches_theshold,aes(theshold,n_matches))+geom_line()

write.csv(manual_checks,"data/full_search/manual_deduplication_test.csv")

dedup_edit<-read.csv("data/full_search/manual_deduplication_test_edit.csv")
head(dedup_edit)
dedup_edit%>%group_by(n)%>%
  summarise(sum_true=sum(true_duplicate==TRUE),sum_false=sum(true_duplicate==FALSE))%>%
  group_by(n)%>%
  summarise(prop_true=sum_true/(sum_true+sum_false))%>%
  ggplot(aes(n,prop_true))+geom_line()

dedup_edit%>%group_by(string_length)%>%
  summarise(sum_true=sum(true_duplicate==TRUE),sum_false=sum(true_duplicate==FALSE))%>%
  group_by(string_length)%>%
  summarise(prop_true=sum_true/(sum_true+sum_false))%>%
  ggplot(aes(string_length,prop_true))+geom_line()

#are additional duplicates that are identified actually duplicates?
#how does this vary by threshold

dedup_edit%>%group_by(n,repeated)%>%
  summarise(sum_true=sum(true_duplicate==TRUE),sum_false=sum(true_duplicate==FALSE))%>%
  ggplot(aes(n,sum_true,colour=repeated))+geom_line()

# in this example, we will use string distance to identify likely duplicates
duplicates_string <- find_duplicates(
  df$title,
  method = "string_osa",
  to_lower = TRUE,
  rm_punctuation = TRUE,
  threshold = 7
)


manual_checks <- review_duplicates(df$title, duplicates_string)

write.csv(manual_checks,"data/full_search_2/manual_deduplication_check.csv")

#plot string length against likelihood to be a duplicate

edited_duplicates<-read.csv("data/full_search/manual_deduplication_check.csv")

head(edited_duplicates)

ggplot(edited_duplicates,aes(length.of.string,score))+geom_point()+geom_smooth(se=F,method = "glm", method.args = list(family = "binomial"))
# it seems pretty obvious that shorter titles result in more erroneous duplicates
#being identified - I will use a string length of 50 as a threshold to check titles

#override incorrectly identified duplicates
unique_matches<-edited_duplicates%>%filter(duplicate=="no")%>% select(matches)%>%summarise(match_unique=unique(matches))


new_duplicates<-override_duplicates(duplicates_string,c(16,260,313,727,1024,1238,3036,3585,3796,3968,4119,4696,
                                        4853,4871,5017,7170,9132,9541))

#remove true duplicates
results <- extract_unique_references(df, new_duplicates)

#set random order for dataframe
set.seed(42)
rows <- sample(nrow(results))
results_random <- results[rows, ]

#create loop to export references in 10% samples
breaks<-data.frame(from=c(1,1337+1,(1337*2)+1,(1337*3)+1,(1337*4)+1,(1337*5)+1,(1337*6)+1,(1337*7)+1,(1337*8)+1,(1337*9)+1),
                   to=c(1337,1337*2,1337*3,1337*4,1337*5,1337*6,1337*7,1337*8,1337*9,13377))
for (n in 1:10){
rand_selection<-results_random[breaks$from[n]:breaks$to[n],]
rand_bib<-write_refs(rand_selection, format="bib",
             file = paste("data/full_search/deduplicated_synthesisr_",n,".bib",sep=""))
}

#export results with duplicates removed
results_ris<-write_refs(results, format="ris",file = "data/full_search/deduplicated_synthesisr.ris")
results_bib<-write_refs(results, format="bib",file = "data/full_search/deduplicated_synthesisr.bib")

