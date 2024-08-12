#Plot the 0fold/4fold dievrsity ratios as boxplots
#The ratios wereobatined with the DiversityRatio.R script
#Open and structure the data for the supegren region, the genome and the 3 supergene arrangements

RatiosSupergene <-read.table(file="Inversion_Ratios_100kb.txt",sep="", header=T)
RatiosSupergene$GenomeRegion<-rep("Supergene",nrow(RatiosSupergene))
RatiosSupergene <-na.exclude(RatiosSupergene) #Remove windows with no ratio
#
RatiosSupergeneHom1<-read.table(file="Separate_InversionHom1_Ratios_100kb.txt",sep="", header=T)
RatiosSupergeneHom1$GenomeRegion<-rep("Chrysippus Homozygote",nrow(RatiosSupergeneHom1))
RatiosSupergeneHom1 <-na.exclude(RatiosSupergeneHom1)#Remove windows with no ratio
colnames(RatiosSupergeneHom1) <-c("scaffold","start","end","sites","zero_four_ratio_Inv","GenomeRegion")

RatiosSupergeneHom2<-read.table(file="Separate_InversionHom2_Ratios_100kb.txt",sep="", header=T)
RatiosSupergeneHom2$GenomeRegion<-rep("Mediterranean Homozygote",nrow(RatiosSupergeneHom2))
RatiosSupergeneHom2 <-na.exclude(RatiosSupergeneHom2)#Remove windows with no ratio
colnames(RatiosSupergeneHom2) <-c("scaffold","start","end","sites","zero_four_ratio_Inv","GenomeRegion")
RatiosSupergeneHet<-read.table(file="Separate_InversionHet_Ratios_100kb.txt",sep="", header=T)
RatiosSupergeneHet$GenomeRegion<-rep("Heterozygote",nrow(RatiosSupergeneHet))
RatiosSupergeneHet <-na.exclude(RatiosSupergeneHet)#Remove windows with no ratio
colnames(RatiosSupergeneHet) <-c("scaffold","start","end","sites","zero_four_ratio_Inv","GenomeRegion")

RatiosGenome <-read.table(file="data/ZeroFour_Ratios_100kb.txt",sep="", header=T)
RatiosGenome$GenomeRegion<-rep("Genome",nrow(RatiosGenome))
RatiosGenome <-na.exclude(RatiosGenome)
#Remove Inversion Windows
RatiosGenome <-RatiosGenome[-c(700:710),]
#RepeatWindows <-which(RatiosGenome$sites == RatiosInversion$sites)

Ratios<-rbind(RatiosGenome,RatiosSupergene)
Ratios3Supergene<-rbind(RatiosSupergeneHom1,RatiosSupergeneHom2,RatiosSupergeneHet)
Ratios3Supergene$LogRatioInver<-log(Ratios3Supergene$zero_four_ratio_Inv) #Log transform the ratios
Ratios$LogRatioMed<-log(Ratios$zero_four_ratio_Med)
Ratios$LogRatioAfr<-log(Ratios$zero_four_ratio_Afr)
#Ratios<-Ratios[c(-850,-927),]
#T-test comparing the ratio between genomes of different regions, between genome and supergene of Mediterranean samples, and supergene of different regions
t.test(Ratios[Ratios$GenomeRegion == "Genome",5],Ratios[Ratios$GenomeRegion == "Genome",6])
t.test(Ratios[Ratios$GenomeRegion == "Supergene",5],Ratios[Ratios$GenomeRegion == "Genome",5])
t.test(Ratios[Ratios$GenomeRegion == "Supergene",5],Ratios[Ratios$GenomeRegion == "Supergene",6])

#Compare the ratios for 3 supergene arrangements
model1 <-  lm(zero_four_ratio_Inv ~ GenomeRegion,data=Ratios3Supergene)
summary(model1)
plot(model1)
#Create a modified data
library(reshape2) 

data_mod <- melt(Ratios, id.vars='GenomeRegion',  
                  measure.vars=c('LogRatioAfr','LogRatioMed'))

data_modInvers <- melt(Ratios3Supergene, id.vars='GenomeRegion',  
                 measure.vars=c('LogRatioInver'))

data_mod_NotLog <- melt(Ratios, id.vars='GenomeRegion',  
                        measure.vars=c('zero_four_ratio_Med','zero_four_ratio_Afr')) 

library(ggplot2)
library(paletteer)
#Plot the results
levels(data_mod$variable) <- c('African','Mediterranean')

data_modInvers$variable <-c("Supergene Genotypes")
data_modInvers$GenomeRegion <- factor(data_modInvers$GenomeRegion , levels=c("Chrysippus Homozygote", "Heterozygote", "Mediterranean Homozygote"))

Asterisks <- data.frame(variable = c("Mediterranean", "African","Supergene Genotypes"),
                       value = c(0.75, 0.75,0))
Asterisks2 <- data.frame(variable = c(1.4, 1.6,2.95),
                        value = c(1.1, 1.6,-1))

r <- 0.15
t <- seq(0, 180, by = 1) * pi / 180
x <- r * cos(t)
y <- r*3 * sin(t)
arc.df <- data.frame(variable = x, value = y)
r <- 0.1
t <- seq(0, 180, by = 1) * pi / 180
x <- r * cos(t)
y <- r*3 * sin(t)
arc.df2 <- data.frame(variable = x, value = y)
line<-data.frame(variable=seq(from=1,to=2,by=0.1),value=2)

ggplot(data_mod, aes(x=variable, y=value)) +
  #scale_y_continuous(breaks= seq(0,260,20)) +
  theme_bw() +
  theme(legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 14),
        legend.key.size =unit(0.5, 'cm')) +
    geom_boxplot(width=0.4,size=0.9,aes(colour=GenomeRegion))+
    geom_boxplot(data= data_modInvers,width=0.4,size=0.9,aes(x=variable,y=value,colour=GenomeRegion))+
    ylab(expression(Log~(pi[0]/pi[4])))+
  xlab("")+
  geom_text(data = Asterisks, label = c("ns","ns","*")) +
  geom_text(data = Asterisks2, label = c("ns","ns","ns")) +
  geom_line(data = arc.df, aes(variable+1, value+0.15), lty = 2) +
  geom_line(data = arc.df, aes(variable+2, value+0.15), lty = 2) +
  geom_line(data = arc.df, aes(variable+3, value-0.55), lty = 2) +
  geom_line(data = arc.df2, aes(variable+2.94, value-1.4), lty = 2) +
  geom_line(data = line, aes(variable+0.1, value-0.5)) +
  geom_line(data = line, aes(variable-0.10, value-1)) +
  scale_colour_paletteer_d("vangogh::Chaise")+
  scale_y_continuous(limits= c(-5,2)) #For not logged 0-2.5, for Logged -5 2
  ggsave("results/100kbRatiosResultsLoged3Invers_v1.pdf",width=15,height=15,units="cm")
  #scale_fill_manual(col=c("#5b942f","#2f6f94"))
