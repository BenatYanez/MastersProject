#Join the HIGH Impact results outputed as one file per chromosome into a single file for the entire genome
#Run in the Cluster
#Position<-0
#chromosomes<-c()

#HighImpact<-data.frame()
#for (i in c(2:20,22:31)){
  #Data<-read.csv(paste("/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/HighImpactVariants.chr",i,".tsv",sep=""),sep="\t")
  #HighImpact<-rbind(HighImpact,Data)
  #TempData<-Data[which(Data$pop == "Mediterranian" ),c(1,2)]
  #TempPosition<-TempData$pos+max(Position)
  #Tempchr<-TempData$chr
  #Position<-c(Position,TempPosition)
 # chromosomes<-c(chromosomes,Tempchr)
#}
#Position<-Position[-1]
#PositionInfor<-data.frame(chromosomes,Position)
#write.table(PositionInfor,file="/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/PositionOfHighImpactVariants.tsv",quote=F,sep="\t",row.names = F)
#write.table(HighImpact,file="/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/HighImpactVariantsGenome.tsv",quote=F,sep="\t",row.names = F)

HighImpact100kb<-data.frame()
for (i in c(2:20,22:31)){
  Data<-read.csv(paste("HighImpactRatio.100kbWindowsSupergene.chr",i,".tsv",sep=""),sep="\t",header=F)
  HighImpact100kb<-rbind(HighImpact100kb,Data)
print(i)
}
write.table(HighImpact100kb,file="/data/martin/genomics/analyses/Danaus_popgen/Benat_project/data/SnpEff/HighImpact100kbWindowsSupergene.tsv",quote=F,sep="\t",row.names=F)
