#Make bed files for 0fold and 4fold sites within the genome and supergene.
#Open the bed file containing all the 0 fold and 4 fold sites in Danaus chrysippus, Obtianed by X code
zero_sites <-read.table(file="data/output.0Dsites.bed")
four_sites <-read.table(file="data/output.4Dsites.bed")
#Only get those sites within the supergene
zero_sites_supergene <-zero_sites[which(zero_sites[,1] == "contig15.1" & zero_sites[,3] <= 6220875 & zero_sites[,3] >= 5322257 | zero_sites[,1] == "contig15.1" & zero_sites[,3] <=7829094  & zero_sites[,3] >= 6251636),]
four_sites_supergene <-four_sites[which(four_sites[,1] == "contig15.1" & four_sites[,3] <= 6220875 & four_sites[,3] >= 5322257 | four_sites[,1] == "contig15.1" & four_sites[,3] <=7829094  & four_sites[,3] >= 6251636),]
#Only get those sites in the rest of the genome
zero_sites_genome <-zero_sites[-which(zero_sites[,1] == "contig15.1" & zero_sites[,3] <= 6220875 & zero_sites[,3] >= 5322257 | zero_sites[,1] == "contig15.1" & zero_sites[,3] <=7829094  & zero_sites[,3] >= 6251636),]
four_sites_genome <-four_sites[-which(four_sites[,1] == "contig15.1" & four_sites[,3] <= 6220875 & four_sites[,3] >= 5322257 | four_sites[,1] == "contig15.1" & four_sites[,3] <=7829094  & four_sites[,3] >= 6251636),]

#Save results
write.table(zero_sites_genome, "output.0Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(four_sites_genome, "output.4Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(zero_sites_supergene, "output.inversion.0Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(four_sites_supergene, "output.inversion.4Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)

   