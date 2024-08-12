#Plot the PCA for the supergene region. The results were obtained with X code
PCA<-read.table(file="data/pca.eu40.vs.dchry2.3.Inversion.csv",sep=",",header=T)
colnames(PCA)[1]<-("Samples")
library(ggplot2)
#Plot the 1st and 2nd proncipal components agianst each other, as they explain the highest variance
ggplot(PCA, aes(x=PC1, y=PC2,col=Samples)) +
  geom_point(size=3) +
  scale_colour_manual(values=c("#3F858C","#3F858C","#3F858C","#707322","#3F858C","#3F858C","#707322","#707322","#707322","#3F858C","#707322","#3F858C","#F2D43D","#707322","#3F858C","#F2D43D","#707322","#F2D43D","#707322","#3F858C","#3F858C","#3F858C","#707322","#3F858C","#707322","#3F858C","#707322","#707322"))+
  theme_bw()+
  theme(legend.position = c(.99, .65),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm'))+
  ylab("PC2 6.2%") +
  xlab("PC1 37.4%")
  #scale_colour_manual(values=c("#FF6978","#3F858C","#B1EDE8","#6D435A","#352D39","#2DE1FC","#2AFC98","#09E85E","#16C172","#214F4B","#30292F","#413F54","#5F5AA2","#355691","#3F4045","#C0C999","#FD96A9","#F62DAE","#B30089","#470063","#343434","#2F3061","#0E34A0","#5F5980","#DFDFDF","#4C6085","#39A0ED","#32322C"))+
ggsave("results/PCA_v1.pdf",width=15,height=15, units="cm")
