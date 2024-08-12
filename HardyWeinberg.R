DiagnosticSites<-read.table(file="data/dchry2.3.chr15.PCRAmplification.3Genotypes.vcf",sep="\t", header=T)
DiagnosticSites<-DiagnosticSites[,-c(1,3,6,7,8,9)]
library(dplyr)
library(tidyr)

SampleNames<-colnames(DiagnosticSites)[-c(1,2,3)]
for (i in 4:33) {
  
  DiagnosticSites<-separate(data=DiagnosticSites,col=i, into = c(paste(colnames(DiagnosticSites)[i],"Genotype",sep=""),NA), 
                            sep = ":", remove = T)
}
j<-1
for (x in (colnames(DiagnosticSites)[c(-1,-2,-3)])) {
  
  DiagnosticSites<-separate(data=DiagnosticSites,col=c(x), into = c(paste(SampleNames[j],"C1",sep=""), paste(SampleNames[j],"C2",sep="")), 
                            sep = "/", remove = FALSE)
  j<-j+1
}

DiagnosticSites<-separate(data=DiagnosticSites,col=c("ALT"), into = c("Allele2", "Allele3"), 
                          sep = ",", remove = TRUE)

for (i in 1:nrow(DiagnosticSites)) {
  DiagnosticSites[i,which(DiagnosticSites[i,] == 1 )] <-DiagnosticSites[i,3]
  DiagnosticSites[i,which(DiagnosticSites[i,] == 0 )]<-DiagnosticSites[i,2]
  DiagnosticSites[i,which(DiagnosticSites[i,] == 2 )]<-DiagnosticSites[i,4]
  
}
#DiagnosticSites_v2<-DiagnosticSites[,c(1,6,7,9,10,12,13,15,16)]
write.table(DiagnosticSites,file="DiagnosticSitesAllSamples.txt")
library(reshape2) 

data_mod <- melt(DiagnosticSites, id.vars='POS',  
                 measure.vars=colnames(DiagnosticSites)[-c(1:4,seq(5,94,by=3))])
data_mod[which(data_mod$value == "."),3]<-NA

colnames(data_mod)<- c("POS","SampleChromosome","Base")
data_mod_SM19SY01 <- melt(DiagnosticSites, id.vars='V2',  
                          measure.vars=c('SM17FV05C1', 'SM17FV05C2',"SM19SY06C1","SM19SY06C2","RV22449C1","RV22449C2","RV12N317C1","RV12N317C2","SM16N01C1","SM16N04C1","RV22450C1","RV22450C2","RV22452C1","RV22452C2"))
data_mod_SM19SY01[which(data_mod_SM19SY01$value == "."),3]<-NA

library(ggplot2)

ggplot(data_mod,aes(x=POS,y=SampleChromosome))+
  geom_point(aes(size=SampleChromosome,alpha=0.5,shape=Base,fill=Base,colour=Base)) +
  theme_bw()+
  xlab("Site")+
  ylab("Sample")+
  scale_colour_manual(values=c("#3F858C","#707322","#F2D43D","#D9814E","blue"))+
  scale_size_manual(values=rep(3,length(unique(data_mod$SampleChromosome)))) +
  guides(size=FALSE, alpha=F)+
  scale_shape_manual(values = c(15,16,17,18,25))
ggsave("results/SNPPCR_AllSamples.png",width=13,height=8)
unique(data_mod$SampleChromosome)
ggplot(data_mod[which(data_mod$SampleChromosome == "RV22449C1" | data_mod$SampleChromosome == "RV22449C2"),],aes(x=POS,y=SampleChromosome))+
  geom_point(aes(size=SampleChromosome,alpha=0.5,shape=Base,fill=Base,colour=Base)) +
  theme_bw()+
  xlab("Site")+
  ylab("Sample")+
  scale_colour_manual(values=c("#3F858C","#707322","#F2D43D","#D9814E","blue"))+
  #xlim(6301891,6301957)+
  xlim(6301957,6301966)+
  scale_size_manual(values=rep(3,length(unique(data_mod$SampleChromosome)))) +
  guides(size=FALSE, alpha=F)+
  scale_shape_manual(values = c(15,16,17,18,25))
data_mod[data_mod$POS<6301957 & data_mod$POS > 6301891,]
ggsave("results/SNPSinPCR_EnhancedRegion_AllSamples.png",width=13,height=8)
########Hardy Weinberg#####
#Two Alleles, Mediterranean African Based on PCA
HomMed<-3+6+2 #Second number is from genotyping PCR
HomAfr<-11+16+16
Het<-12+20+6 #Not entirely sure these numbers are correct
Individuals<-HomMed+HomAfr+Het
PObsFreq<-(2*(HomAfr)+Het)/(2*(HomMed+HomAfr+Het))
QObsFreq<-1-PObsFreq
ExpHomAfr<-PObsFreq^2*Individuals
ExpHomMed<-QObsFreq^2*Individuals
ExpHet<-2*PObsFreq*QObsFreq*Individuals
Expected<-c(ExpHomAfr,ExpHet,ExpHomMed)
Observed<-c(HomAfr,Het,HomMed)
Results<-data.frame(Observed,Expected)
Test<-sum(((Observed-Expected)^2)/Expected) #Chis-qr stat 0.00249

#Three Alleles, Mediterranean, Chrysippus, Orientis-Like  
HomMed<-3+6+2
HomOr<-3+9+11
HomChr<-3+6+1
HetMedOr<-7+17+5
HetMedChr<-5+3+1
HetOrChr<-5+1+4

Individuals<-sum(HomMed+HomOr+HomChr+HetMedOr+HetMedChr+HetOrChr)
MedObsFreq<-(2*(HomMed)+HetMedOr+HetMedChr)/(2*Individuals)
OrObsFreq<-(2*(HomOr)+HetMedOr+HetOrChr)/(2*Individuals)
ChrObsFreq<-(2*(HomChr)+HetMedChr+HetOrChr)/(2*Individuals)
ExpHomMed<-MedObsFreq^2*Individuals
ExpHomOr<-OrObsFreq^2*Individuals
ExpHomChr<-ChrObsFreq^2*Individuals
ExpHetMedOr<-2*MedObsFreq*OrObsFreq*Individuals
ExpHetMedChr<-2*MedObsFreq*ChrObsFreq*Individuals
ExpHetOrChr<-2*OrObsFreq*ChrObsFreq*Individuals
Expected<-c(ExpHomMed,ExpHomOr,ExpHomChr,ExpHetMedOr,ExpHetMedChr,ExpHetOrChr)
Observed<-c(HomMed,HomOr,HomChr,HetMedOr,HetMedChr,HetOrChr)
Results<-data.frame(Observed,Expected)
Test<-chisq.test(Results) #0.8789 Not significant
Test<-sum((((HomMed-ExpHomMed)^2)/ExpHomMed)+(((HomOr-ExpHomOr)^2)/ExpHomOr)+(((HomChr-ExpHomChr)^2)/ExpHomChr)+(((HetMedOr-ExpHetMedOr)^2)/ExpHetMedOr)+(((HetMedChr-ExpHetMedChr)^2)/ExpHetMedChr)+(((HetOrChr-ExpHetOrChr)^2)/ExpHetOrChr))
Test<-sum(((Observed-Expected)^2)/Expected) #13.769
#2 df
library(HardyWeinberg)
HWChisq(y)
plot(x=DiagnosticSites$V2,y=as.numeric(DiagnosticSites$SM17FV05C1))
