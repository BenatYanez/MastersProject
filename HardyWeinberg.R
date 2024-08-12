
#Using the SNPs from the PCR amplified region, identify diagnostic sites to differentiate the Chysippus,Mediterranean and Orientis allele using genome data
#Open file of diagnostic SNPs obtained with the X code
DiagnosticSites<-read.table(file="data/dchry2.3.chr15.PCRAmplification.3Genotypes.vcf",sep="\t", header=T)
DiagnosticSites<-DiagnosticSites[,-c(1,3,6,7,8,9)]
library(dplyr)
library(tidyr)
#Just obtain the names of the samples from the columns in the VCF
SampleNames<-colnames(DiagnosticSites)[-c(1,2,3)]
#Remove all data except for genotypes from each cell
for (i in 4:33) {
  
  DiagnosticSites<-separate(data=DiagnosticSites,col=i, into = c(paste(colnames(DiagnosticSites)[i],"Genotype",sep=""),NA), 
                            sep = ":", remove = T)
}
#Separate 0/1 genotypes into two columns , one column with 0 and another with 1
j<-1
for (x in (colnames(DiagnosticSites)[c(-1,-2,-3)])) {
  
  DiagnosticSites<-separate(data=DiagnosticSites,col=c(x), into = c(paste(SampleNames[j],"C1",sep=""), paste(SampleNames[j],"C2",sep="")), 
                            sep = "/", remove = FALSE)
  j<-j+1
}
#Separate the alternate allele into 2 columns, so that if the site is tri-allelic the Alelle3 columns has that base 
DiagnosticSites<-separate(data=DiagnosticSites,col=c("ALT"), into = c("Allele2", "Allele3"), 
                          sep = ",", remove = TRUE)

#Substitute the 0,1,2 with bases. A genotype of 0 represents the Reference allele, 1 is Allele2 and 2 is Allele3
for (i in 1:nrow(DiagnosticSites)) {
  DiagnosticSites[i,which(DiagnosticSites[i,] == 1 )] <-DiagnosticSites[i,3]
  DiagnosticSites[i,which(DiagnosticSites[i,] == 0 )]<-DiagnosticSites[i,2]
  DiagnosticSites[i,which(DiagnosticSites[i,] == 2 )]<-DiagnosticSites[i,4]
  
}
#Save the Diagnostic sites
write.table(DiagnosticSites,file="results/DiagnosticSitesAllSamples.txt")
#Preparet the dataframe to plot
library(reshape2) 
data_mod <- melt(DiagnosticSites, id.vars='POS',  
                 measure.vars=colnames(DiagnosticSites)[-c(1:4,seq(5,94,by=3))])
data_mod[which(data_mod$value == "."),3]<-NA

colnames(data_mod)<- c("POS","SampleChromosome","Base")
#data_mod_SM19SY01 <- melt(DiagnosticSites, id.vars='V2',  
                #          measure.vars=c('SM17FV05C1', 'SM17FV05C2',"SM19SY06C1","SM19SY06C2","RV22449C1","RV22449C2","RV12N317C1","RV12N317C2","SM16N01C1","SM16N04C1","RV22450C1","RV22450C2","RV22452C1","RV22452C2"))
#data_mod_SM19SY01[which(data_mod_SM19SY01$value == "."),3]<-NA

library(ggplot2)
#Plot all the sites
ggplot(data_mod,aes(x=POS,y=SampleChromosome))+
  geom_point(aes(size=SampleChromosome,alpha=0.5,shape=Base,fill=Base,colour=Base)) +
  theme_bw()+
  xlab("Site")+
  ylab("Sample")+
  scale_colour_manual(values=c("#3F858C","#707322","#F2D43D","#D9814E","blue"))+
  scale_size_manual(values=rep(3,length(unique(data_mod$SampleChromosome)))) +
  guides(size=FALSE, alpha=F)+
  scale_shape_manual(values = c(15,16,17,18,25))
ggsave("results/SNPPCR_AllSamples.pdf",width=15,height=9,units="cm")

#Plot only a particular region, which can differentiate the 3 alleles
#ggplot(data_mod[which(data_mod$SampleChromosome == "RV22449C1" | data_mod$SampleChromosome == "RV22449C2"),],aes(x=POS,y=SampleChromosome))+
ggplot(data_mod,aes(x=POS,y=SampleChromosome))+
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
ggsave("results/SNPSinPCR_EnhancedRegion_AllSamples.pdf",width=15,height=9,units="cm")
########Hardy Weinberg#####
#Using the allele frequencies from from the genomes and using Sanger sequencing calculate deviations from Hardy Weinberg
#Two Alleles, Mediterranean African, Looking at the Supergene as a whole
HomMed<-3+6+2 #First number is from genomes, 2nd and 3rd from Sanger sequences
HomAfr<-11+16+16
Het<-12+20+6 
Individuals<-HomMed+HomAfr+Het
PObsFreq<-(2*(HomAfr)+Het)/(2*(HomMed+HomAfr+Het))
QObsFreq<-1-PObsFreq
ExpHomAfr<-PObsFreq^2*Individuals
ExpHomMed<-QObsFreq^2*Individuals
ExpHet<-2*PObsFreq*QObsFreq*Individuals
Expected<-c(ExpHomAfr,ExpHet,ExpHomMed)
Observed<-c(HomAfr,Het,HomMed)
Results<-data.frame(Observed,Expected)
Test<-sum(((Observed-Expected)^2)/Expected) #Chis-qr stat 0.3335784

#Three Alleles, Mediterranean, Chrysippus, Orientis-Like. Just looking at the B locus  
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
Test<-sum(((Observed-Expected)^2)/Expected) #13.769
#2 df
