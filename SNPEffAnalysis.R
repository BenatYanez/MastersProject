#Obtain  the ratio of HIGH impact variants to Callable sites obtained by the how_high.py script
#Plot the ratios for sueprgene, genome of the two regions and for the 3 supergene arrangements and comapre them with t-tests

#Open data
HighImpactWindows<-read.csv("data/HighImpact100kbWindowsGenome.tsv",sep="\t")
SupergeneHighImpactWindows<-read.csv("data/HighImpact100kbWindowsSupergene.tsv",sep="\t")
colnames(HighImpactWindows)<-c("Chromosome","Start","End","HighImpact","Callable","Location")
colnames(SupergeneHighImpactWindows)<-c("Chromosome","Start","End","HighImpact","Callable","Arrangement")
#Obtain ratios of HIGH impact to CAllable sites
HighImpactWindows$Ratio<-HighImpactWindows$HighImpact/HighImpactWindows$Callable
SupergeneHighImpactWindows$Ratio<-SupergeneHighImpactWindows$HighImpact/SupergeneHighImpactWindows$Callable
#Separate the windows that occur within the supergene
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


SupergeneHighImpactWindows$GenomeRegion<-SupergeneHighImpactWindows$Region
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Genome",]) #Between the 3 datastes the only difference appears to be between the Mediterranean and African samples for the entire genome
#T-test comparing the ratio between genomes of different regions, between genome and supergene of Mediterranean samples, and supergene of different regions
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Genome",])
t.test(Ratio ~ Region,data=HighImpactWindows[HighImpactWindows$Location == "Mediterranian",])
t.test(Ratio ~ Region,data=HighImpactWindows[HighImpactWindows$Location == "African",])
t.test(Ratio ~ Location,data=HighImpactWindows[HighImpactWindows$Region == "Supergene",])
#Compare the ratios for 3 supergene arrangements
model1 <-  lm(Ratio ~ Arrangement,data=SupergeneHighImpactWindows )
summary(model1)
plot(model1)
library(paletteer)
library(ggplot2)
#Plot boxplot of the ratios
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
  ggsave(paste("results/HighImpactVariantsBoxplotTruncated.pdf",sep=""),width=15,height=15,units="cm")



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
