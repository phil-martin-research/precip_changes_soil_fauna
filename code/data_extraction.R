#script for data extraction of a few papers

library(metaDigitise)
library(shinyDigitise)


metaDigitise("data/data_extraction/Ashton_2019/")

Ashton_2019_data<-metaDigitise(dir="data/data_extraction/Ashton_2019/")
write.csv(Ashton_2019_data,file = "data/data_extraction/Ashton_2019/Ashton_2019_figures.csv")


#extract data for Zhou_2023
Zhou_data<-shinyDigitise("C:/Users/phil.martin/Desktop/Zhou_2023")

