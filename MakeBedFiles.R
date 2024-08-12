zero_sites <-read.table(file="data/output.0Dsites.bed")
four_sites <-read.table(file="data/output.4Dsites.bed")

zero_sites_inversion <-zero_sites[which(zero_sites[,1] == "contig15.1" & zero_sites[,3] <= 6220875 & zero_sites[,3] >= 5322257 | zero_sites[,1] == "contig15.1" & zero_sites[,3] <=7829094  & zero_sites[,3] >= 6251636),]
four_sites_inversion <-four_sites[which(four_sites[,1] == "contig15.1" & four_sites[,3] <= 6220875 & four_sites[,3] >= 5322257 | four_sites[,1] == "contig15.1" & four_sites[,3] <=7829094  & four_sites[,3] >= 6251636),]

zero_sites_genome <-zero_sites[-which(zero_sites[,1] == "contig15.1" & zero_sites[,3] <= 6220875 & zero_sites[,3] >= 5322257 | zero_sites[,1] == "contig15.1" & zero_sites[,3] <=7829094  & zero_sites[,3] >= 6251636),]
four_sites_genome <-four_sites[-which(four_sites[,1] == "contig15.1" & four_sites[,3] <= 6220875 & four_sites[,3] >= 5322257 | four_sites[,1] == "contig15.1" & four_sites[,3] <=7829094  & four_sites[,3] >= 6251636),]


write.table(zero_sites_genome, "output.0Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(four_sites_genome, "output.4Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(zero_sites_inversion, "output.inversion.0Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(four_sites_inversion, "output.inversion.4Dsites.bed",sep="\t",col.names=F,row.names=F,quote=F)

   