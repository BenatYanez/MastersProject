
library(dplyr)
library(tidyr)

HighImpact<-read.csv("data/HighImpactVariantsGenome.tsv",sep="\t")

Position<-read.csv("data/PositionOfHighImpactVariants.tsv",sep="\t")

HighImpact<-separate(data=HighImpact,col="n_high_alleles", into = c("HighEffect1", "HighEffect2"), 
           sep = ",", remove = FALSE)



MedHighImpact<-HighImpact[which(HighImpact$pop == "Mediterranian" ),c(1,2,3,4,6,7)]
AfrHighImpact<-HighImpact[which(HighImpact$pop == "African" ),c(1,2,3,4,6,7)]
MedHighImpact$HighEffect1<-as.numeric(MedHighImpact$HighEffect1)
AfrHighImpact$HighEffect1<-as.numeric(AfrHighImpact$HighEffect1)
MedHighImpact$HighEffect2<-as.numeric(MedHighImpact$HighEffect2)
AfrHighImpact$HighEffect2<-as.numeric(AfrHighImpact$HighEffect2)
MedHighImpact<-MedHighImpact %>% mutate(HighEffect2 = ifelse(is.na(HighEffect2), 0, HighEffect2))
AfrHighImpact<-AfrHighImpact %>% mutate(HighEffect2 = ifelse(is.na(HighEffect2), 0, HighEffect2))

MedHighImpact$HighFreq<-((MedHighImpact$HighEffect1 + MedHighImpact$HighEffect2)/MedHighImpact$n_alleles)
AfrHighImpact$HighFreq<-((AfrHighImpact$HighEffect1 + AfrHighImpact$HighEffect2)/AfrHighImpact$n_alleles)

#AfrHighImpact$HighFreq<-AfrHighImpact$HighEffect1/AfrHighImpact$n_alleles #Only include the 1st allele 
#MedHighImpact$pos<-Position$Position
#AfrHighImpact$pos<-Position$Position

MedResultsWindows100kb<-data.frame()
AfrResultsWindows100kb<-data.frame()
for (i in unique(HighImpact$chr)) { 
  MedWindows100kb<-data.frame()
  AfrWindows100kb<-data.frame()
  Windows<-seq(from=1, to=max(HighImpact[which(HighImpact$chr == i),2]),by=100000)
  MedChrmData<-MedHighImpact[which(MedHighImpact$chr == i),]
  AfrChrmData<-AfrHighImpact[which(AfrHighImpact$chr == i),]
  
  for (x in 1:c(length(Windows)-1)) {
    
  MedWindows100kb[x,1]<-Windows[x]
  MedWindows100kb[x,2]<-Windows[x+1]-1
  MedWindows100kb[x,3]<-mean(MedChrmData[MedChrmData$pos < Windows[x+1] & MedChrmData$pos >= Windows[x],7],na.rm=T)
  MedWindows100kb[x,4]<-nrow(MedChrmData[MedChrmData$pos < Windows[x+1] & MedChrmData$pos >= Windows[x] & MedChrmData$HighFreq != 0,])
  MedWindows100kb[x,5]<-i

  
  AfrWindows100kb[x,1]<-Windows[x]
 AfrWindows100kb[x,2]<-Windows[x+1]
  AfrWindows100kb[x,3]<-mean(AfrChrmData[AfrChrmData$pos < Windows[x+1] & AfrChrmData$pos >= Windows[x],7],na.rm=T)
  AfrWindows100kb[x,4]<-nrow(AfrChrmData[AfrChrmData$pos < Windows[x+1] & AfrChrmData$pos >= Windows[x],])
  AfrWindows100kb[x,5]<-i
  
  }
  MedResultsWindows100kb<-rbind(MedResultsWindows100kb,MedWindows100kb)
  AfrResultsWindows100kb<-rbind(AfrResultsWindows100kb,AfrWindows100kb)
  
}
colnames(MedResultsWindows100kb)<-c("Start","End","HighFreq","Sites","Chromosome")
colnames(AfrResultsWindows100kb)<-c("Start","End","HighFreq","Sites","Chromosome")
MedResultsWindows100kb$Midpoint<-seq(from=0.5,by=1,length.out=nrow(MedResultsWindows100kb))
AfrResultsWindows100kb$Midpoint<-seq(from=0.5,by=1,length.out=nrow(AfrResultsWindows100kb))

tickpos<-c()
for (i in unique(MedResultsWindows100kb$Chromosome)) {
  if (i == "chr02") { tickpos<-max(MedResultsWindows100kb[MedResultsWindows100kb$Chromosome == i,6])/2}
  else {tickpos<-c(tickpos,((max(MedResultsWindows100kb[MedResultsWindows100kb$Chromosome == i,6])-min(MedResultsWindows100kb[MedResultsWindows100kb$Chromosome == i,6]))/2)+min(MedResultsWindows100kb[MedResultsWindows100kb$Chromosome == i,6]))}
}
Label<-c(2:20,22:31)
MedResultsWindows100kb$Test<-1/(MedResultsWindows100kb$HighFreq/MedResultsWindows100kb$Sites)
library(ggplot2)
library(paletteer)

  ggplot(MedResultsWindows100kb,aes(x=Midpoint,y=HighFreq))+
  geom_point(aes(col=Chromosome,size=Sites)) +
  # geom_point(data=InversionWindowsAfr,shape=15,alpha=0.5,aes(x=Midpoint,y=HighFreq))+
  #geom_point(data=InversionWindows,shape=17,alpha=0.5,aes(x=Midpoint,y=HighFreq))+
  
  #geom_point(data=AfrResultsWindows100kb,shape=15,aes(x=Midpoint,y=HighFreq,alpha=0.5))+
  scale_x_continuous(name="Chromosome", breaks=tickpos,labels=Label) +
  scale_color_manual(values=c(rep(c("#F2D43D","#731A12"),14),"#F2D43D"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=14))
ggsave("InversionsHighlighted100kbWindowsHighFreq_v2.png",width=16,height= 8.55)
MedResultsWindows100kb$Region<-"Mediterranean"
AfrResultsWindows100kb$Region<-"Africa"

InversionWindows<-MedResultsWindows100kb[MedResultsWindows100kb$Chromosome == "chr15" & MedResultsWindows100kb$Start >= 5322257 & MedResultsWindows100kb$End < 6220875 | MedResultsWindows100kb$Chromosome == "chr15" & MedResultsWindows100kb$Start >= 6251636 & MedResultsWindows100kb$End < 7829094,  ]
InversionWindowsAfr<-AfrResultsWindows100kb[AfrResultsWindows100kb$Chromosome == "chr15" & AfrResultsWindows100kb$Start >= 5322257 & AfrResultsWindows100kb$End < 6220875 | AfrResultsWindows100kb$Chromosome == "chr15" & AfrResultsWindows100kb$Start >= 6251636 & AfrResultsWindows100kb$End < 7829094,  ]
InversionWindows$Region<-"Mediterranean"
InversionWindowsAfr$Region<-"Africa"
StatData<-cbind(InversionWindows,InversionWindowsAfr)
NonInversionData<-MedResultsWindows100kb[-(MedResultsWindows100kb$Chromosome == "chr15" & MedResultsWindows100kb$Start >= 5322257 & MedResultsWindows100kb$End < 6220875 | MedResultsWindows100kb$Chromosome == "chr15" & MedResultsWindows100kb$Start >= 6251636 & MedResultsWindows100kb$End < 7829094),  ]
NonInversionData$Area<-"Genome"
InversionWindows$Area<-"Inversion"
StatData2<-rbind(NonInversionData,InversionWindows)
#Tests whether the Windows inside the "Inversion" have a greater frequency than those outside
t.test(x=InversionWindows$HighFreq,y=NonInversionData$HighFreq,alternative="greater") #Significantly higher in Inversion
#Test whetehr the inversion in the Mediterranean population has a higher frequency of HIgh impact variants than the African samples
t.test(x=StatData[,3],y=StatData[,10],paired=T,alternative="greater",data=StatData) #Significantly higher in Mediterranean
#Maybe something non parametric, weight

ggplot(InversionWindows,aes(x=Midpoint,y=HighFreq))+
  geom_point() +
  #geom_point(data=AfrResultsWindows100kb,aes(x=Midpoint,y=HighFreq,col="green"))+
  scale_x_continuous(name="Chromosome", breaks=tickpos,labels=Label) +
  scale_color_manual(values=c(rep(c("blue","red"),14),"blue"))+
  theme_bw()+
  theme(legend.position = "none",axis.text.x = element_text(size=14))


#Windows of High Impact based on Thomas how_high script
HighImpactWindows<-read.csv("data/HighImpact100kbWindowsGenome.tsv",sep="\t")
SupergeneHighImpactWindows<-read.csv("data/HighImpact100kbWindowsSupergene.tsv",sep="\t")
colnames(HighImpactWindows)<-c("Chromosome","Start","End","HighImpact","Callable","Location")
colnames(SupergeneHighImpactWindows)<-c("Chromosome","Start","End","HighImpact","Callable","Arrangement")
HighImpactWindows$Ratio<-HighImpactWindows$HighImpact/HighImpactWindows$Callable
SupergeneHighImpactWindows$Ratio<-SupergeneHighImpactWindows$HighImpact/SupergeneHighImpactWindows$Callable
SupergeneHighImpactWindows <-SupergeneHighImpactWindows[SupergeneHighImpactWindows$Chromosome == "chr15" & SupergeneHighImpactWindows$Start >= 5322257 & SupergeneHighImpactWindows$End < 6220875 | SupergeneHighImpactWindows$Chromosome == "chr15" & SupergeneHighImpactWindows$Start >= 6251636 & SupergeneHighImpactWindows$End < 7829094,]
SupergeneHighImpactWindows$Region<-"Supergene"
SupergeneHighImpactWindows[SupergeneHighImpactWindows$Arrangement == "Het",6]<-"Heterozygote"
SupergeneHighImpactWindows[SupergeneHighImpactWindows$Arrangement == "Hom1",6]<-"Chrysippus Homozygote"
SupergeneHighImpactWindows[SupergeneHighImpactWindows$Arrangement == "Hom2",6]<-"Mediterranean Homozygote"
SupergeneHighImpactWindows$Arrangement <- factor(SupergeneHighImpactWindows$Arrangement , levels=c("Chrysippus Homozygote", "Heterozygote", "Mediterranean Homozygote"))

HighImpactWindows$Region<-"Genome"
HighImpactWindows$GenomeRegion<-HighImpactWindows$Region

HighImpactWindows[HighImpactWindows$Chromosome == "chr15" & HighImpactWindows$Start >= 5322257 & HighImpactWindows$End < 6220875 | HighImpactWindows$Chromosome == "chr15" & HighImpactWindows$Start >= 6251636 & HighImpactWindows$End < 7829094,8]<-"Supergene"

HighImpactWindows$GenomeRegion<-HighImpactWindows$Region

#Remove windows with less than 1000 callable sites #This is the same as removing windows where having just two HIGh impact alleles would result in a value as high as the 3rd Quantile
HighImpactWindowsRL1000<-HighImpactWindows[HighImpactWindows$Callable>= 1000,]
hist(HighImpactWindowsRL1000$Callable,breaks=1000,xlim=c(0,50000))
summary(HighImpactWindows$Ratio)
#Remove windows with less than 2000 callable sites #This is the same as removing windows where having just four HIGh impact alleles would result in a value as high as the 3rd Quantile
HighImpactWindowsRL2000<-HighImpactWindows[HighImpactWindows$Callable>= 2000,]
hist(HighImpactWindowsRL2000$Callable,breaks=1000,xlim=c(0,50000))
summary(HighImpactWindowsRL2000$Ratio)
#Remove windows with less than 10000
HighImpactWindowsRL10000<-HighImpactWindows[HighImpactWindows$Callable>= 10000,]
hist(HighImpactWindowsRL10000$Callable,breaks=1000,xlim=c(0,50000))
plot(y=HighImpactWindowsRL10000$Ratio,x=HighImpactWindowsRL10000$Callable)
plot(y=HighImpactWindows$HighImpact,x=HighImpactWindows$Callable)
points(y=HighImpactWindows[which(HighImpactWindows$Region == "Supergene"),4],x=HighImpactWindows[which(HighImpactWindows$Region == "Supergene"),5],pch=15,col="red")
summary(HighImpactWindowsRL10000$Ratio)
SupergeneHighImpactWindowsRL10000<-SupergeneHighImpactWindows[SupergeneHighImpactWindows$Callable>= 10000,]
model1
SupergeneHighImpactWindows$GenomeRegion<-SupergeneHighImpactWindows$Region
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Genome",]) #Between the 3 datastes the only difference appears to be between the Mediterranean and African samples for the entire genome
#T tests
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Genome",])
t.test(Ratio ~ Region,data=HighImpactWindows[HighImpactWindows$Location == "Mediterranian",])
t.test(Ratio ~ Region,data=HighImpactWindows[HighImpactWindows$Location == "African",])
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Supergene",])
model2 <-  lm(Ratio ~ Arrangement,data=SupergeneHighImpactWindows )
summary(model2)
plot(model2)

#This might be because the MEditerranean sampels have overall a lower NE than African ones so there is an accumulation of deleterious variants. This could indicate that there is not enough power for the supergene region.
#Removing windows of less than 10000 does not have an impact in the different haplotypes for the supergene

ggplot(HighImpactWindows, aes(x=Location, y=Ratio)) +
  #scale_y_continuous(breaks= seq(0,260,20)) +
  theme_bw() +
  theme(legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = ), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
  geom_boxplot(width=0.4,size=0.9,aes(colour=GenomeRegion))+
  geom_boxplot(data= SupergeneHighImpactWindows,width=0.4,size=0.9,aes(x=Region,y=Ratio,colour=Arrangement)) +
  ylab("Number of High Impact variants / Callable Variants")+
  xlab("")+
  scale_colour_paletteer_d("vangogh::Chaise")+
  scale_y_continuous(limits= c(0,0.005))
  ggsave(paste("results/TestHighImpactVariantsBoxplotTruncated.pdf",sep=""),width=15,height=15,units="cm")

#Remove windows with less than 1000 callable sites
hist(HighImpactWindows$Callable,breaks=1000,xlim=c(0,50000))
summary(HighImpactWindows$Ratio)
plot(y=HighImpactWindows$Ratio,x=HighImpactWindows$Callable)
#Try to cut the lower tail possibly, try several thresholds,




InversionRegion<- HighImpactWindows[HighImpactWindows$Chromosome == "chr15" & HighImpactWindows$Start >= 5322257 & HighImpactWindows$End < 6220875 | HighImpactWindows$Chromosome == "chr15" & HighImpactWindows$Start >= 6251636 & HighImpactWindows$End < 7829094,]



MedHighImpactWindows<-HighImpactWindows[HighImpactWindows$Location == "Mediterranian",]
MedHighImpactWindows$Midpoint<-seq(from=0.5,by=1,length.out=nrow(MedHighImpactWindows))
AfrHighImpactWindows<-HighImpactWindows[HighImpactWindows$Location == "African",]
AfrHighImpactWindows$Midpoint<-seq(from=0.5,by=1,length.out=nrow(AfrHighImpactWindows))

model<-glm(Ratio ~ Location*Region,data=HighImpactWindows)
anova(model)
summary(model)
plot(model)
library(dplyr)
library(tidyr)

Data_Mode<-data.frame(Start=unique(SupergeneHighImpactWindows$Start),Het=SupergeneHighImpactWindows[seq(1,nrow(SupergeneHighImpactWindows),by=3),7],Hom1=SupergeneHighImpactWindows[seq(2,nrow(SupergeneHighImpactWindows),by=3),7],Hom2=SupergeneHighImpactWindows[seq(3,nrow(SupergeneHighImpactWindows),by=3),7])
t.test(Data_Mode$Hom1,Data_Mode$Hom2,paired=T)
Data_Mode2<-data.frame(Start=unique(InversionRegion$Start),African=InversionRegion[seq(1,nrow(InversionRegion),by=2),7],Meditarranian=InversionRegion[seq(2,nrow(InversionRegion),by=2),7])
t.test(Data_Mode2$African,Data_Mode2$Meditarranian,paired=T)

data_mod<- melt(SupergeneHighImpactWindows, id.vars='Start',  
                measure.vars=c("Arrangement"))

library(paletteer)

ggplot(HighImpactWindows, aes(x=Location, y=Ratio)) +
  #scale_y_continuous(breaks= seq(0,260,20)) +
  theme_bw() +
  theme(legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 22),
        legend.key.size =unit(1, 'cm')) +
  geom_boxplot(width=0.4,size=0.9,aes(colour=GenomeRegion))+
  geom_boxplot(data= SupergeneHighImpactWindows,width=0.4,size=0.9,aes(x=GenomeRegion,y=Ratio,colour=Arrangement)) +
  ylab("Number of High Impact variants / Callable Variants")+
  xlab("")+
  scale_colour_paletteer_d("vangogh::Chaise") + 
  scale_y_continuous(limits= c(0,0.005))

ggsave(paste("results/TestHighImpactVariantsBoxplotAll.png",sep=""),width=8,height=8)
ggsave(paste("results/TestHighImpactVariantsBoxplotAll.pdf",sep=""),width=8,height=8)


tickpos<-c()
for (i in unique(MedHighImpactWindows$V1)) {
  if (i == "chr02") { tickpos<-max(MedHighImpactWindows[MedHighImpactWindows$V1 == i,8])/2}
  else {tickpos<-c(tickpos,((max(MedHighImpactWindows[MedHighImpactWindows$V1 == i,8])-min(MedHighImpactWindows[MedHighImpactWindows$V1 == i,8]))/2)+min(MedHighImpactWindows[MedHighImpactWindows$V1 == i,8]))}
}
Label<-c(2:20,22:31)
library(ggplot2)
ggplot(MedHighImpactWindows,aes(x=Midpoint,y=Ratio))+
  geom_point(aes(col=V1)) +
  # geom_point(data=InversionWindowsAfr,shape=15,alpha=0.5,aes(x=Midpoint,y=HighFreq))+
  #geom_point(data=InversionWindows,shape=17,alpha=0.5,aes(x=Midpoint,y=HighFreq))+
  geom_point(data=AfrHighImpactWindows,shape=15,aes(x=Midpoint,y=Ratio,alpha=0.5))+
  scale_x_continuous(name="Chromosome", breaks=tickpos,labels=Label) +
  scale_color_manual(values=c(rep(c("#F2D43D","#731A12"),14),"#F2D43D"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=14))
