library(jpeg)
library("exiftoolr")
library("tidyverse")

######
TemperatureData<-data.frame()
v<-1
CenterX  <- c(651,652,653)
CenterY<- c(480,481,482)
#Each day I measured temperature had different background temperature, upper limits and lower limits. 
#So each batch establishes which individuals I ran, which data they were measured and the abckgrounda nd limit temperatures
#To run the code craete the vectors for one Batch and then run the For loop bellow, then do the next batch an run the loop again 
#June 3th Background 22ºC Limits 18.4ºC-38ºC
Date<- c("03_06_2024")
Batch<- c(9,19,12,5)
BackgroundReal<-22
Upper<-38
Lower<-18.4
#June 3th Background 22ºC Limits 10.9ºC-30.3ºC
Batch<- c(26,35,14,5,1,28,6)
Upper<-30.3
Lower<-10.9
#June 4th Background 22.1ºC Limits 10.6ºC-30ºC
Date<- c("04_06_2024")
Batch<- c(18,8)
BackgroundReal<-22.1
Upper<-30
Lower<-10.6
#June 4th Background 24.2ºC Limits 10.6ºC-30ºC
Batch<- c(33,21,27,11,2,22,7,17,31,16,30,34,13,5,31,28,35,19,12,15,14,30,36,26,1,23,4,9,10,3)
BackgroundReal<-24.2
#June 5th Background 24.7ºC Limits 10.8-30.4ºC
Date<- c("05_06_2024")
Batch  <- c(17)
BackgroundReal<-24.7
Upper<-30.4
Lower<-10.8
#June 5th Background 25.5ºC Limits 10.8-30.4ºC
Batch  <- c(6,32,21,18,2,16,22,8,11,27,33,14,1,12,9,5,19,13,4,35,10,36,20,26,30,23,31)
BackgroundReal<-25.5
#June 5th Background 23ºC Limits 10.8-30.4ºC
Batch<- c("33_2","11_2","27_2","2_2","22_2","3_2","18_2","32_2","6_2","17_2","8_2","21_2","16_2")
BackgroundReal<-23
#June 6th Background 25ºC Limits 10.2ºC-30.6ºC
Date<- c("06_06_2024")
Batch<- c(26,10,14,35,23,13,9,31,36,5,12,19,2,33,8,18,21,16,17,3,11,20,32,1,22,4,27)
BackgroundReal<-25
Upper<-30.6
Lower<-10.2
#######
#The data consists of folders for the 4 days, with another folder for each individual measured that day
#Each folder then has all the thermal images for one run
Range<-Upper-Lower
for (j in Batch) { 
  my_dir <- path.expand(paste("data/ThermalImages/",Date,"/",j,"/",sep="")) #Determine where the photos are for each batch
  file_names<- list.files(path=my_dir, pattern = "*.JPG")
  creation_dates <- paste0(my_dir,file_names)    %>% #The next few lines obtain the time the picture was taken from its metadata
    lapply(function(x) exif_read(x)$FileModifyDate) %>% #the timepoint for the first picture is time 0, then for the next add the time between pictures
    unlist
  ModifyDates<- sub("01:00","",creation_dates)
  ModifyDates<- sub("[+]","",ModifyDates)
  ModifyDates<- sapply(strsplit(ModifyDates, " "), function(x) x[2])         
  Times<- sapply(strsplit(ModifyDates, ":"), function(n) as.numeric(n) %*% c(3600, 60, 1))
  n<- 0
  Temp<-data.frame()

  for (x in file_names) { 
    image<- readJPEG(paste(my_dir,x,sep="")) #Open a picture as a raster
    CenterImage<-mean(image[CenterY,CenterX,]) #Obtain the pixel value of the center
    Background<-Lower+Range*median(image[162:801,332:971,]) #Obtain the tempeature of the Background from its pixel value
    Correction<-BackgroundReal-Background #Correct the temperature based on the difference between infered and real background temperature
    Temp[x,1]<-Lower+Range*CenterImage+ Correction #Calculate temperature of the center
    Temp[x,2]  <- Times[1+n]-Times[1]
    n<-n+1
    Temp[x,3]<-j
    Temp[x,4]<-Date
    Temp[x,5]<-v
  }
 v<-v+1
    TemperatureData<-rbind(TemperatureData,Temp)
}

TemperatureData$Sample<-gsub("_2","",TemperatureData$Sample)
#Save results
colnames(TemperatureData)<-(c("Temperature","Time","Sample","Date","Run"))
write.csv(TemperatureData,"results/Results.csv")
####### 
#Clugii data, Different morph
#June 6th Background 25ºC Limits 10.2ºC-30.6ºC
Date<- c("06_06_2024")
Batch<- c("1","2","3","4","2_2","3_2","4_2")
BackgroundReal<-25
Upper<-30.6
Lower<-10.2
#######
Range<-Upper-Lower
ClugiiData<-data.frame()

for (j in Batch) { 
  my_dir <- path.expand(paste("data/ThermalImages/",Date,"/","Clugii",j,"/",sep=""))
  file_names<- list.files(path=my_dir, pattern = "*.JPG")
  creation_dates <- paste0(my_dir,file_names)    %>% 
    lapply(function(x) exif_read(x)$FileModifyDate) %>% 
    unlist
  ModifyDates<- sub("01:00","",creation_dates)
  ModifyDates<- sub("[+]","",ModifyDates)
  ModifyDates<- sapply(strsplit(ModifyDates, " "), function(x) x[2])         
  Times<- sapply(strsplit(ModifyDates, ":"), function(n) as.numeric(n) %*% c(3600, 60, 1))
  n<- 0
  Temp<-data.frame()
  
  for (x in file_names) { 
    image<- readJPEG(paste(my_dir,x,sep=""))
    CenterImage<-mean(image[CenterY,CenterX,])
    Background<-Lower+Range*median(image[162:801,332:971,])
    Correction<-BackgroundReal-Background
    Temp[x,1]<-Lower+Range*CenterImage+ Correction
    Temp[x,2]  <- Times[1+n]-Times[1]
    n<-n+1
    Temp[x,3]<-j
    Temp[x,4]<-Date
    Temp[x,5]<-v
  }
  v<-v+1
  ClugiiData<-rbind(ClugiiData,Temp)
}

colnames(ClugiiData)<-(c("Temperature","Time","Sample","Date","Run"))

###############
