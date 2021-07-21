#Cleaning working space
rm(list=ls()) 
#Creator: 20160817 FPL
#Modifications: 20170408 CG
today<-format(Sys.Date(),"%Y%m%d")
#--------------------------------------------
#From Linux
# get working directory
#setwd("D:/Dropbox/FL+WB/PELA/data")
# Cristóbal
setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones BIMAC/Plasticidad pela/analisis/datos")
getwd()     #verifico directorio
#----------------------------------------------

#-----------------------------------------------
#Libraries
#-----------------------------------------------
library(plotrix);library(ggplot2);library(nlme);
library(car);library(multcomp)
#-----------------------------------------------

#------------------------------------------------
#reading my data
#------------------------------------------------
pela<-read.csv("data.pela.csv", header=TRUE, 
               sep=",",dec=".", strip.white=TRUE)
names(pela)
summary(pela)
str(pela)

#transformo variables mal asignadas a numericas
pela$stage=as.factor(pela$stage)
pela$oxygen=as.numeric(pela$oxygen)
pela$saturation=as.numeric(pela$saturation)
pela$mo2=as.numeric(pela$mo2)
pela$temperature=as.factor(pela$temperature)
# defino kpa promedio para cada tratamiento
pela$kpa.mean<-ifelse(pela$kpa==2.07 | pela$kpa==2.36 | pela$kpa==2.65, 2.36, 
                      ifelse(pela$kpa==4.14 | pela$kpa==4.71 | pela$kpa==5.30, 4.71,
                             ifelse(pela$kpa==8.29 | pela$kpa==9.43 | pela$kpa==10.60, 9.44,
                                    ifelse(pela$kpa==12.43 | pela$kpa==14.14 | pela$kpa==15.89, 14.15,
                                           ifelse(pela$kpa==18.65 | pela$kpa==21.21 | pela$kpa==23.84, 21.23, NA)))))

table(pela$stage, pela$kpa.mean)
# ordeno niveles etapa
pela$stage<- factor(pela$stage, levels = c("Embryo","Megalopae","Juvenile","Adult"))
# Estimo media y sd
pela.mean<-aggregate(cbind(mo2.3,mo2.4)~stage+kpa.mean+temperature,data=pela, mean)
pela.sd<-aggregate(cbind(mo2.3,mo2.4)~stage+kpa.mean+temperature,data=pela, sd)
pela1<-merge(pela.mean, pela.sd, by=c("stage","kpa.mean","temperature"))
names(pela1)[c(4:7)]<-c("mo2.3", "mo2.4", "mo2.3.sd", "mo2.4.sd")
#---------------------------------------------------
# Consuno de oxygeno por etapa y temperatura
#---------------------------------------------------
{
png(filename="Figure S1.png",width=18,height=18,units="cm",res=600)
par(mfrow=c(4,3),oma=c(3,1,1,1),tcl=-0.3, family="Arial", mar=c(2,4.5,0,0)) #mar fija margenes alrededor de figuras
#--------------------------------------
#1.0.Embryos
#--------------------------------------
with(subset(pela1,temperature==6 & stage=="Embryo"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),las=1,ylab="", xaxt="n"))
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
#mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=1)
text(10,39,c("Eggs 6°C acclimated"), cex=0.9)

with(subset(pela1,temperature==12 & stage=="Embryo"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), xaxt="n",yaxt="n",las=1, ylab=""))
text(10,39,c("Eggs 12°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)

with(subset(pela1,temperature==18 & stage=="Embryo"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), xaxt="n",yaxt="n",las=1, ylab=""))
text(10,39,c("Eggs 18°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=-2, outer=T)
#--------------------------------------
#2.0.Megalopae
#--------------------------------------
with(subset(pela1,temperature==6 & stage=="Megalopae"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),las=1,xaxt="n",ylab=""))
#mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=1)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
text(12,39,c("Megalopae 6°C acclimated"), cex=0.9)

with(subset(pela1,temperature==12 & stage=="Megalopae"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),xaxt="n" ,yaxt="n",las=1, ylab=""))
text(12,39,c("Megalopae 12°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
#text(20,c(0.09),c("12°C acclimated"))
with(subset(pela1,temperature==18 & stage=="Megalopae"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), xaxt="n",yaxt="n",las=1, ylab=""))
text(12,39,c("Megalopae 18°C acclimated"), cex=0.9)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
#--------------------------------------
#3.0.Juveniles
#--------------------------------------
with(subset(pela1,temperature==6 & stage=="Juvenile"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),las=1,xaxt="n",ylab=""))
text(10,39,c("Juveniles 6°C acclimated"), cex=0.9)
#mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=1)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)

with(subset(pela1,temperature==12 & stage=="Juvenile"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), xaxt="n",yaxt="n",las=1, ylab=""))
text(10,39,c("Juveniles 12°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)

with(subset(pela1,temperature==18 & stage=="Juvenile"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),xaxt="n" ,yaxt="n",las=1, ylab=""))
text(10,39,c("Juveniles 18°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
axis(1,at = seq(5,20,5), labels = F, tick = TRUE, lwd.ticks = 0.6)
#--------------------------------------
#4.0.Juveniles
#--------------------------------------
with(subset(pela1,temperature==6 & stage=="Adult"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40),las=1,ylab=""))
abline(a=NULL,b=NULL,v=9.16, lty=2)
text(10,39,c("Adults 6°C acclimated"), cex=0.9)
#mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=1)

with(subset(pela1,temperature==12 & stage=="Adult"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), yaxt="n",las=1, ylab=""))
abline(a=NULL,b=NULL,v=15.46, lty=2)
text(10,39,c("Adults 12°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)

with(subset(pela1,temperature==18 & stage=="Adult"),plotCI(kpa.mean,mo2.3, mo2.3.sd, ylim=c(0,40), yaxt="n",las=1, ylab=""))
abline(a=NULL,b=NULL,v=15.66, lty=2)
text(10,39,c("Adults 18°C acclimated"), cex=0.9)
axis(2,at = seq(0,40,10), labels = F, tick = TRUE, lwd.ticks = 0.6)
mtext("Oxygen tension (kPa) ",side=1,  line=1, outer=T, at=0.55)

dev.off()
}
#-----------------------------------------------------------------------------------------------------------------
#  2da Revisión
#------------------------------------------------------------------------------------------------------------------
#Reviewer 2: I would suggest to depict for each of the four life stages an additional figure S2 (A-D) 
#            presenting the effect of temperature (x-axis) on oxygen consumption rates (y-axis) 
#            for each oxygen tension in one graph.
pela1$temperature<-as.numeric(as.character(pela1$temperature))
{
png(filename="Figure S3..png",width=14,height=12,units="cm",res=600)
par(mfrow=c(2,2),oma=c(1,2,1,0),tcl=-0.3, family="Arial", mar=c(2,2,0,1)) #mar fija margenes alrededor de figuras
  
# embryos
with(subset(pela1,kpa.mean==2.36 & stage=="Embryo"),plotCI(temperature,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="white", ylim=c(0,40), xlim=c(6,24)))
with(subset(pela1,kpa.mean==4.71 & stage=="Embryo"),plotCI(temperature+1,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="gray", add=T))
with(subset(pela1,kpa.mean==9.44 & stage=="Embryo"),plotCI(temperature+2,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="black", add=T))
with(subset(pela1,kpa.mean==14.15 & stage=="Embryo"),plotCI(temperature+3,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="white", add=T))
with(subset(pela1,kpa.mean==21.23 & stage=="Embryo"),plotCI(temperature+4,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="gray", add=T))
axis(1,at=c(6,10),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=c(12,16),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=c(18,22),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
text(6,38, labels = "Eggs",     adj = c( 0, 1 ), col = "black" )
# leyenda
legend("topright", legend=c("2.36 KPa", "4.71 kPa","9.44 kPa","14.15 kPa","21.23 kPa"),pch=c(21,21,21,22,22),
       pt.bg=c("white", "gray","black","white","gray"), cex=0.7)
# Megalopae
with(subset(pela1,kpa.mean==2.36 & stage=="Megalopae"),plotCI(temperature,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="white", ylim=c(0,40),xlim=c(6,24)))
with(subset(pela1,kpa.mean==4.71 & stage=="Megalopae"),plotCI(temperature+1,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="gray", add=T))
with(subset(pela1,kpa.mean==9.44 & stage=="Megalopae"),plotCI(temperature+2,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="black", add=T))
with(subset(pela1,kpa.mean==14.15 & stage=="Megalopae"),plotCI(temperature+3,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="white", add=T))
with(subset(pela1,kpa.mean==21.23 & stage=="Megalopae"),plotCI(temperature+4,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="gray", add=T))
axis(1,at=c(6,10),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=c(12,16),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=c(18,22),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
text(6,38, labels = "Megalopae",     adj = c( 0, 1 ), col = "black" )
# Juvenile
with(subset(pela1,kpa.mean==2.36 & stage=="Juvenile"),plotCI(temperature,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="white", ylim=c(0,40),xlim=c(6,24)))
with(subset(pela1,kpa.mean==4.71 & stage=="Juvenile"),plotCI(temperature+1,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="gray", add=T))
with(subset(pela1,kpa.mean==9.44 & stage=="Juvenile"),plotCI(temperature+2,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="black", add=T))
with(subset(pela1,kpa.mean==14.15 & stage=="Juvenile"),plotCI(temperature+3,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="white", add=T))
with(subset(pela1,kpa.mean==21.23 & stage=="Juvenile"),plotCI(temperature+4,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="gray", add=T))
text(6,38, labels = "Juveniles",     adj = c( 0, 1 ), col = "black" )
# etiquetas eje x
# identifica estado caligus
axis(1,at=c(6,10),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=8,col="black",line=0,labels="6",tick=T)
axis(1,at=c(12,16),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=14,col="black",line=0,labels="12",tick=T)
axis(1,at=c(18,22),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=20,col="black",line=0,labels="18",tick=T)
# Adult
with(subset(pela1,kpa.mean==2.36 & stage=="Adult"),plotCI(temperature,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="white", ylim=c(0,40),xlim=c(6,24)))
with(subset(pela1,kpa.mean==4.71 & stage=="Adult"),plotCI(temperature+1,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="gray", add=T))
with(subset(pela1,kpa.mean==9.44 & stage=="Adult"),plotCI(temperature+2,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=21,pt="black", add=T))
with(subset(pela1,kpa.mean==14.15 & stage=="Adult"),plotCI(temperature+3,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="white", add=T))
with(subset(pela1,kpa.mean==21.23 & stage=="Adult"),plotCI(temperature+4,mo2.3, mo2.3.sd,las=1,ylab="", xaxt="n",pch=22,pt="gray", add=T))
text(6,38, labels = "Adults",     adj = c( 0, 1 ), col = "black" )
# etiquetas eje x
axis(1,at=c(6,10),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=8,col="black",line=0,labels="6",tick=T)
axis(1,at=c(12,16),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=14,col="black",line=0,labels="12",tick=T)
axis(1,at=c(18,22),col="black",line=0,tick=T,labels=rep("",2),lwd=2,lwd.ticks=0)
axis(1,at=20,col="black",line=0,labels="18",tick=T)
# etiqueta eje y
mtext(expression(paste(dot(M)[O[2]], " (",mu,"mol ", O[2], " ",h^-1, " ", g^-1, ")")),side=2,  line=0.5, outer=T)
mtext("Temperature (°C)",side=1, line=0, outer=T)
dev.off()
}
