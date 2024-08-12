contigs<-read.table("data/contigs.txt",sep="")
chromEnds<-c(12593762,12257276,16603682,9939236,11368998,12307431,17804611,13759957,7589886,4330961,7253754,3510316,6377628,7525859,8495143,2159359,7108895,8413355,5438580,4304363,5109367,4406282,3598968,8903490,15527882,3639899,5430338,5133822,1926300,2486427,7985327,5410907,4431358,16112269,5031374,8326359,12027074,12415182,11492671,13367733,11449205)
  #Want to create 100 "windows" that are 2.5Mb at random Make sure the point is created bellow the end of each chromosome
for (i in 1:nrow(contigs)) {
  WindowStart <- 0 
  WindowEnd <- WindowStart + 2500000
 bedfile <- data.frame("V1"= contigs[i,1],
                "V2"= WindowStart,
                "V3"= WindowEnd)
 x<-1
 
  repeat{ 
    x<-x+1
    WindowStart<- WindowEnd+1
    WindowEnd <- WindowStart + 2500000
    bedfile[x,1]<- contigs[i,1]
    bedfile[x,2]<- WindowStart
    bedfile[x,3]<-WindowEnd
  if (WindowStart >= chromEnds[[i]]-2500000) 
    (break)
  }
 outtext <-paste("results/WindowBeds/Windows of 2.5Mb",contigs[i,1],".bed",sep="")
 write.table(bedfile, file = outtext,
             quote= F,col.names = F)
 }
  
sample (cont)
i <-2
#Create bed files for each chromosome of 2.5Mb long for later analysis of the frequency of high impact variants with SNPEFF
chromEnds<-c(1,12549439,16112267,13360988,12027094,12415243,11492700,13367733,11449226,12593783,12257276,9939308,11369012,12307432,17804635,13761402,11921163,10764171,1390589,10656133,0,8413376,9414270,8005422,8903523,15527911,9070343,7060234,10474926,5410915,4431405)
AutosomeWindows <- data.frame()
#Want to create 100 "windows" that are 2.5Mb at random Make sure the point is created bellow the end of each chromosome
               for (i in c(2:20,22:31)) {
                 WindowStart <- 0
                 WindowEnd <- WindowStart + 2500000
                 y<- i
                 if(i <10) (
                   y<-paste("0",i,sep="")
                 )
                 bedfile <- data.frame("V1"= paste("chr",y,sep=""),
                                       "V2"= WindowStart,
                                       "V3"= WindowEnd)
                 x<-1
                 
                 repeat{
                   x<-x+1
                   WindowStart<- WindowEnd+1
                   if (WindowStart >= chromEnds[[i]]-2500000)
                     (break)
                   WindowEnd <- WindowStart + 2500000

                   bedfile[x,1]<- paste("chr",y,sep="")
                   bedfile[x,2]<- WindowStart
                   bedfile[x,3]<-WindowEnd
                 
                 }
                 AutosomeWindows <- rbind(AutosomeWindows,bedfile)
                 outtext <-paste("results/WindowBeds/Windows_2.5Mb_chr",i,".bed",sep="")
                 write.table(bedfile, file = outtext,
                             quote= F,col.names = F)
                 write.table(AutosomeWindows,file="AutosomeWindows_2.5Mb.bed",quote = F,col.names = F)
               }
