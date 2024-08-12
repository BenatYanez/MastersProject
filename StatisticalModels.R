library(ggplot2)
#Run statistical analysis of the temperature data obtained by ImageTemperature.R script
#Load Data
TemperatureData<-read.csv("results/Results_Temperature.csv",sep=";")
TemperatureData$Run<-as.factor(TemperatureData$Run)


#Add the Wing Morph 
SampleInfo<-read.csv("data/Spain_samples.csv",sep=";")
for (i in 1:nrow(TemperatureData)) {
  TemperatureData[i,7]<-SampleInfo[which(SampleInfo$ï..ID.Number. == TemperatureData[i,4]),3]
  TemperatureData[i,8]<-SampleInfo[which(SampleInfo$ï..ID.Number. == TemperatureData[i,4]),5]
  TemperatureData[i,9]<-SampleInfo[which(SampleInfo$ï..ID.Number. == TemperatureData[i,4]),8]
  TemperatureData[i,10]<-SampleInfo[which(SampleInfo$ï..ID.Number. == TemperatureData[i,4]),9]
  
}
colnames(TemperatureData)<-c("PictureID","Temperature","Time","Sample","Date","Run","WingMorph","AbdomenMorph","WingLength","WingHeight")
TemperatureData$FactorTime<-as.factor(TemperatureData$Time)

#How many runs for each Morph?
Morph<-c()
for(i in 1:111) {
  Morph[i]<- unique(TemperatureData[TemperatureData$Run == i, 7])
}
sum(Morph == "Dark")
sum(Morph == "Pale")

#Cut data so that it starts when butterfly thorax reaches 21º
ModifiedTemperatureData <-data.frame()
for (i in unique(TemperatureData$Run)) {
  if (min(TemperatureData[which(TemperatureData$Run == i),2])>21) {next}
  wheretemp1<- max(which(TemperatureData$Run == i & TemperatureData$Temperature <= 21))
  temp1<- TemperatureData[wheretemp1,2]
  vector<- which(TemperatureData$Run == i & c(TemperatureData$Temperature >= 21 ))
  wheretemp2<- min(vector[vector >= wheretemp1])
  temp2<- TemperatureData[wheretemp2,2]
  time1 <-TemperatureData[which(TemperatureData$Run == i & TemperatureData$Temperature == temp1),3]
  time2<-TemperatureData[which(TemperatureData$Run == i & TemperatureData$Temperature == temp2),3]
  slope<-(time2-time1)/(temp2-temp1)
  time0<- (21/slope)+time1
  times<- TemperatureData[which(TemperatureData$Run == i & TemperatureData$Time >= time2),3]-time0
  temps<- TemperatureData[which(TemperatureData$Run == i & TemperatureData$Time >= time2),2]
  Morph<- unique(TemperatureData[which(TemperatureData$Run == i), 7])
  ID<- unique(TemperatureData[which(TemperatureData$Run == i), 4])
  Date<-unique(TemperatureData[which(TemperatureData$Run == i), 5])
  WingLength<-unique(TemperatureData[which(TemperatureData$Run == i), 9])
  WingHeight<-unique(TemperatureData[which(TemperatureData$Run == i), 10])
  Temporary<-data.frame(c(0,times),c(21,temps),rep(i,length(temps)+1),rep(Morph,length(temps)+1),rep(ID,length(temps)+1),rep(Date,length(temps)+1),rep(WingLength,length(temps)+1),rep(WingHeight,length(temps)+1))
  ModifiedTemperatureData<-rbind(ModifiedTemperatureData,Temporary)
}

colnames(ModifiedTemperatureData)<-c("Time","Temperature","Run","Morph","Sample","Date","WingLength","WingHeight")
ModifiedTemperatureData$FactorTime<-as.factor(ModifiedTemperatureData$Time)
#Plot Raw data
ggplot(TemperatureData,aes(x=Time,y=Temperature,col=WingMorph,group=Run)) +
  geom_point()+
  geom_line(lty=1,alpha=0.5) +
  theme_bw() +
  labs(col = "Wing Morph") +
  theme(legend.position = c(.99, .50),
        legend.justification = c("right"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size = 17), 
        legend.title = element_text(size = 20),
        legend.key.size =unit(1, 'cm')) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ylab("Temperature (ºC)")+
  xlab("Time (s)")+
  scale_color_manual(values=c("#731A12","#F2D43D"))

ggsave("results/RawDataTemp_v2.png",width=12,height=7)

#Scale so that the quadratic is not in a broadly different scale
#Do not mean center otherwise the intercept will not be correct
TemperatureData_Scaled<-TemperatureData
TemperatureData_Scaled$Time<-scale(TemperatureData$Time,center=F)
TemperatureData_Scaled$WingLength<-scale(TemperatureData$WingLength,center=F)
TemperatureData_Scaled$WingHeight<-scale(TemperatureData$WingHeight,center=F)

library(glmmTMB)
library(ggeffects)
#Model with lowest AIC score
#ar1 covariance structure should account for autocorrelation due to repeated measurements
#Quadratic Model
TemperatureData_Scaled$WingMorph<-as.factor(TemperatureData_Scaled$WingMorph)
LowestAIC_QuadModel<-glmmTMB(Temperature ~Date+Time*WingMorph+I(Time^2)*WingMorph+(1+Time+I(Time^2)|Run)+ ar1(FactorTime + 0 | Run),data=TemperatureData_Scaled,family="gaussian")
#Obtain predictions based on the model
LowestAIC_QuadModel_Predict<-ggpredict(LowestAIC_QuadModel,terms=c("Time[all]","WingMorph","Date[04_06_2024]"),interval="confidence",type="random")
#Plot predictions
ggplot(LowestAIC_QuadModel_Predict,aes(x=x,y=predicted,colour=group))+
  geom_line()+
  theme_bw() +
  ylab("Temperature (ºC)")+
  xlab("Time (s)")+
  labs(colour="Wing Morph",fill="") +
  scale_x_continuous(breaks=c(0,100/135.9927,200/135.9927,300/135.9927,400/135.9927), #135.9927 is the standard deviation of the original dataset used to scale down time
                     labels = ~ round((135.9927*(.x)))) +
  theme(legend.position = c(.30, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
  scale_color_manual(values=c("#731A12","#F2D43D")) +
  scale_fill_manual(values=c("#731A12","#F2D43D")) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), alpha = 0.1)
ggsave("results/FullDataQuadModelEffects_GGPredict.pdf",width=15,height=9,units="cm")

summary(LowestAIC_QuadModel) #5963 AIC
ConfidenceIntervals<- confint(WingSize_LogModel)
ConfidenceIntervals[1:12,] #Obtain confidence intervals for the estimates

#Use a model where I log Transform Time
WingSize_LogModel<-glmmTMB(Temperature ~Date*log1p(Time)+log1p(Time)*scale(WingHeight)+log1p(Time)*scale(WingLength)+log1p(Time)*WingMorph+(1|Run)+ar1(FactorTime + 0 | Run),data=TemperatureData)
WingSize_LogModel_Predict<-ggpredict(WingSize_LogModel,terms=c("Time[all]","WingMorph","Date[04_06_2024]"),interval="confidence",type="random")
par(mfrow = c(1, 3))
ggplot(WingSize_LogModel_Predict,aes(x=x,y=predicted,colour=group))+
  geom_line()+
  theme_bw() +
  scale_x_continuous(trans = "log1p")+
  ylab("Temperature (ºC)")+
  xlab("Time (s)")+
  labs(colour="Wing Morph",fill="") +
  theme(legend.position = c(.30, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
  scale_color_manual(values=c("#F2D43D","#731A12")) +
  scale_fill_manual(values=c("#F2D43D","#731A12")) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), alpha = 0.1)
ggsave("results/FullDataLogModelEffects_GGPredict.pdf",width=15,height=9,units="cm")

summary(WingSize_LogModel)#AIC 9453

ConfidenceIntervals<- confint(WingSize_LogModel)
ConfidenceIntervals[1:15,]


#Run the same Models but on the 21º Modified Data
ModifiedTemperatureData_Scaled<-ModifiedTemperatureData
ModifiedTemperatureData_Scaled$Time<-scale(ModifiedTemperatureData$Time,center=F)
ModifiedTemperatureData_Scaled$WingHeight<-scale(ModifiedTemperatureData$WingHeight,center=F)
ModifiedTemperatureData_Scaled$WingLength<-scale(ModifiedTemperatureData$WingLength,center=F)


LowestAIC_Quad21Model<-glmmTMB(Temperature ~Date*Time+Time*Morph+I(Time^2)*Morph+(1+Time+I(Time^2)|Run)+ ar1(FactorTime + 0 | Run),data=ModifiedTemperatureData_Scaled)
summary(LowestAIC_Quad21Model) #AIC 3601
ConfidenceIntervals<- confint(LowestAIC_Quad21Model)
ConfidenceIntervals[1:12,]
LowestAIC_Quad21Model_Predict<-ggpredict(LowestAIC_Quad21Model,terms=c("Time[all]","Morph","Date[04_06_2024]"),interval="confidence",type="random")
ggplot(LowestAIC_Quad21Model_Predict,aes(x=x,y=predicted,colour=group))+
  geom_line()+
  theme_bw() +
  ylab("Temperature (ºC)")+
  xlab("Time (s)")+
  labs(colour="Wing Morph",fill="") +
  scale_x_continuous(breaks=c(0,100/135.9927,200/135.9927,300/135.9927,400/135.9927), #135.9927 is the standard deviation of the original dataset used to scale down time
                     labels = ~ round((135.9927*(.x)))) +
  theme(legend.position = c(.30, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
  scale_color_manual(values=c("#F2D43D","#731A12")) +
  scale_fill_manual(values=c("#F2D43D","#731A12")) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), alpha = 0.1)
ggsave("results/21ºCDataQuadModelEffectsGGPredict.pdf",width=15,height=9,units="cm")
#Run Log transformed model
ModData_LowestAIC_LogModel2<-glmmTMB(Temperature ~Date*log1p(Time)+log1p(Time)*scale(WingHeight)+log1p(Time)*scale(WingLength)+log1p(Time)*Morph+(1|Run)+ar1(FactorTime + 0 | Run),data=ModifiedTemperatureData)
ModData_LowestAIC_LogModel2_Predict<-ggpredict(ModData_LowestAIC_LogModel2,terms=c("Time[all]","Morph","Date[04_06_2024]"),interval="confidence",type="random")
ggplot(ModData_LowestAIC_LogModel2_Predict,aes(x=x,y=predicted,colour=group))+
  geom_line()+
  theme_bw() +
  scale_x_continuous(trans = "log1p")+
  ylab("Temperature (ºC)")+
  xlab("Time (s)")+
  labs(colour="Wing Morph",fill="") +
  theme(legend.position = c(.30, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
  scale_color_manual(values=c("#F2D43D","#731A12")) +
  scale_fill_manual(values=c("#F2D43D","#731A12")) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), alpha = 0.1)
ggsave("results/21ºCDataLogModelEffectsGGEffect.pdf",width=15,height=9,unit="cm")

summary(ModData_LowestAIC_LogModel2)#4384 AIC

ConfidenceIntervals<- confint(ModData_LowestAIC_LogModel2)
ConfidenceIntervals[1:15,]

#Add impact of date
#Model with log1p



summary(WingSize_LogModel)
#Check assumptions in DHARMA
library(DHARMa)
testDispersion(WingSize_Model)
simulationOutput <- simulateResiduals(fittedModel = LowestAIC_QuadModel, plot = F)


plot(simulationOutput)

