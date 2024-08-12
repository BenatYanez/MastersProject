#Carry out a power analysis of the ability to detect devaitions from hardy Weinberg due to Heterozygote advanatge

q<-0.3^0.5 #Light morph frequency BB, based the frequency of a few pale morph
p<-1-q
Genotypes<-c("AA","BA","BB")
#Before selection
Heterozygote<-2*p*q
Homozygoteq<-q^2
Homozygotep<-p^2
#Vector to store the pvalue
pvalue<-numeric(length=10000)
#Different sample size
SampleSizes<-seq(from=20, to=200, by=20)
#Different strengths of selection
S<-seq(from=0, to=1, by=0.05)
#Store the results
PowerResults<-data.frame(matrix(NA,nrow=10,ncol=21))
colnames(PowerResults)<-S
rownames(PowerResults)<-SampleSizes
#Run Power analysis for each sample size and selection strength
for (j in 1:10){
  
  for(x in 1:21){
   #Selection impacts both homozygote sthe same, reducing their frequency
     SelHeterozygote<-2*p*q
    SelHomozygoteq<-q^2*(1-S[x])
    SelHomozygotep<-p^2*(1-S[x])
    
    for (i in 1:10000) {
#Sample both homozygotes and the heterozygotes with probability with no Heterozygote advantage
      Samples<- sample(Genotypes,size=SampleSizes[j],replace=T,prob=c(Homozygotep,Heterozygote,Homozygoteq))
      SampledAB<-length(which(Samples == "BA"))
      SampledBB<-length(which(Samples == "BB"))
      SampledAA<-length(which(Samples == "AA"))
      #Estimate the allele frequencies from the samples
      Estimatep<-(2*SampledAA+SampledAB)/(2*SampleSizes[j])
      Estimateq<- 1-Estimatep
      #Estimate what the genotype frequencies would be under Hardy Weinberg
      ExpectedAA<-(Estimatep^2)*SampleSizes[j]
      ExpectedAB<-Estimatep*Estimateq*2*SampleSizes[j] 
      ExpectedBB<-(Estimateq^2)*SampleSizes[j]
      
      #Simulate getting individuals that have undergone selection for a heterozygote with strength s
      SelSamples<- sample(Genotypes,size=SampleSizes[j],replace=T,prob=c(SelHomozygotep,SelHeterozygote,SelHomozygoteq))
      SelSampledAB<-length(which(SelSamples == "BA"))
      SelSampledBB<-length(which(SelSamples == "BB"))
      SelSampledAA<-length(which(SelSamples == "AA"))
      #Create a dataframe with both
      Expected<-c(ExpectedAA,ExpectedAB,ExpectedBB)
      Observed<-c(SelSampledAA,SelSampledAB,SelSampledBB)
      Results<-data.frame(Expected,Observed)
      
      #Do a chi squared test to see whether these are any different and store the pvalue for the power analysis
      pvalue[i]<-chisq.test(Results)$p.val
    }
    PowerResults[j,x]<-length(which(pvalue <= 0.05 ))/length(pvalue)
    
  }
  
}
#Store and plot results
write.csv(PowerResults, file="results/PowerResults(Sel 1-S).csv")
library(ggplot2)
#Reshape Dataframe
df_reshaped <- data.frame (x = rep(as.numeric(rownames(PowerResults)),21),y=c(PowerResults$`0`,PowerResults$`0.05`,PowerResults$`0.1`,PowerResults$`0.15`,PowerResults$`0.2`,PowerResults$`0.25`,PowerResults$`0.3`,PowerResults$`0.35`,PowerResults$`0.4`,PowerResults$`0.45`,PowerResults$`0.5`,PowerResults$`0.55`,PowerResults$`0.6`,PowerResults$`0.65`,PowerResults$`0.7`,PowerResults$`0.75`,PowerResults$`0.8`,PowerResults$`0.85`,PowerResults$`0.9`,PowerResults$`0.95`,PowerResults$`1`) ,group=as.factor(c(rep(0, 10),rep(0.05, 10),rep(0.1, 10),rep(0.15, 10),rep(0.2, 10),rep(0.25, 10),rep(0.3, 10),rep(0.35, 10),rep(0.4, 10),rep(0.45, 10),rep(0.5, 10),rep(0.55, 10),rep(0.6, 10),rep(0.65, 10),rep(0.7, 10),rep(0.75, 10),rep(0.8, 10),rep(0.85, 10),rep(0.9, 10),rep(0.95, 10),rep(1, 10))))
#Colour palette                
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "#5f7213", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "blue1", 
  "green1", "yellow4",
  "darkorange4", "brown"
)

ggplot(df_reshaped,aes( x=x,y, col=as.factor(group))) +
  geom_line(linewidth=1.2)+
  geom_hline(yintercept=0.8)+
  scale_colour_manual (values=c25)+
  labs(x="Sample Size",y="Power", col="Selection Strength") +
  theme_bw() 
ggsave("results/PowerAnalysis3.pdf",width=15,height=9, units="cm")

