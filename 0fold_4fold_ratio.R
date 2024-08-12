#Calculate the Ratio of 0fold to 4fold diversity in the genome(Autosomes)

#Create a dataframe for the results
Results <- data.frame()
size<-100
#four_diversity_All<-read.table(file=paste("data/diversity/Autosome.w",size,"kb.m1kb.4D",sep=""),sep=",", header=T)
#four_diversity_All$scaffold<-sub("chr0","",four_diversity_All$scaffold)
#four_diversity_All$scaffold<-sub("chr","",four_diversity_All$scaffold)

for (i in c(2:20,22:31)) {
  zero_diversity<-read.table(file=paste("data/diversity/chr",i,".w",size,"kb.m1kb.0D",sep="" ),sep=",", header=T)
  four_diversity<-read.table(file=paste("data/diversity/chr",i,".w",size,"kb.m1kb.4D",sep="" ),sep=",", header=T)
 # four_diversity<-four_diversity_All[which(four_diversity_All$scaffold == i),]
  Ratios<-zero_diversity[,c(1,2,3,5)]
  Ratios$zero_four_ratio_Med<-zero_diversity$pi_Mediterranian/four_diversity$pi_Mediterranian
  Ratios$zero_four_ratio_Afr<-zero_diversity$pi_African/four_diversity$pi_African
  Results<-rbind(Results,na.exclude(Ratios))
}
write.table(Results,file=paste("data/ZeroFour_Ratios_",size,"kb.txt",sep=""))
