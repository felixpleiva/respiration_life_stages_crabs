#---------
#Title
#---------
#Respiratory control and thermal sensibility in benthic 
#life stages of porcelain crab (Petrolisthes laevigatus).
#-----------------------------------------------
#Cleaning working space
rm(list=ls()) 
#Creator: 20160817 FPL
#Modifications: 20170408 CG
today<-format(Sys.Date(),"%Y%m%d")
#--------------------------------------------
#From Linux
# get working directory
# setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
#From Windows
# setwd("D:/felix.leiva@ulagos.cl/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
setwd("C:/Users/Invunche/Dropbox/Radboud University/publicaciones/Publicadas/3. Leiva et al Mar Biol 2018/data")
getwd()     #verifico directorio
#-----------------------------------------------
#Libraries
#-----------------------------------------------
library(nlme);library(car);library(multcomp)
#-----------------------------------------------
#reading my data
#------------------------------------------------
pela<-read.csv("data.pela.csv", header=TRUE, 
               sep=",",dec=".", strip.white=TRUE)
names(pela)
summary(pela)
str(pela)
#------------------------------------------------

#---------------------------------------------------------------------------
#Units for each explicative variable
#---------------------------------------------------------------------------
#dry_weight_g: individual dry weight,expressed in grams
#mo2: metabolic rate, expressed as mg O2 por hour per individual
#mo2.1: metabolic rate, expressed as mg O2 per hour per gram dry weight
#mo2.2: metabolic rate, expressed as mmol O2 per hour per gram dry weight
#mo2.3: metabolic rate, expressed as umol O2 per hour per gram dry weight
#mo2.4: metabolic rate, expressed as umol O2 por hour per individual
#--------------------------------------------------------
#Transformation of variables wrongly assigned to numeric
#--------------------------------------------------------
pela$stage=as.factor(pela$stage)
pela$oxygen_mg=as.numeric(pela$oxygen_mg)
pela$saturation=as.numeric(pela$saturation)
pela$mo2=as.numeric(pela$mo2)
pela$kpa<-as.numeric(pela$kpa)
pela$temperature=as.factor(pela$temperature)
###########################################
#EXPORTO MEDIAS PARA METAZOA
grouped <- group_by(pela,stage,kpa.mean,temperature)
datos.mean<-summarise(grouped, mean=mean(mo2.4), sd=sd(mo2.4))

round_df <- function(datos.mean, digits) {
  nums <- vapply(datos.mean, is.numeric, FUN.VALUE = logical(1))
  datos.mean[,nums] <- round(datos.mean[,nums], digits = digits)
  
  (datos.mean)
}

datos.met<-round_df(datos.mean, digits=6)

data.pela.metazoa<-merge(datos.mass,datos.met,by=c("stage","kpa.mean","temperature"))
data.pela.metazoa<-rename(data.pela.metazoa,dry_mass_mean=mean.x)
data.pela.metazoa<-rename(data.pela.metazoa,dry_mass_sd=sd.x)
data.pela.metazoa<-rename(data.pela.metazoa,resp_mean=mean.y)
data.pela.metazoa<-rename(data.pela.metazoa,resp_sd=sd.y)
#write.csv(data.pela.metazoa,"6090 Supplementary data of respiration.csv",row.names = FALSE)
###########################################
#-------------------------------------------------------------
#Scaling of Metabolic rates
#--------------------------------------------------------
# subset by stage
#--------------------------------------------------------
datos<-pela
unique(pela$stage)
embrio<-subset(datos,datos$stage=="Embryo")
adult<-subset(datos,datos$stage=="Adult")
juvenil<-subset(datos,datos$stage=="Juvenile")
megalopa<-subset(datos,datos$stage=="Megalopae")
Embryo01<-c(1:588)
Embryo01[which(pela$stage!="Embryo")]<-0 #without embryos
Embryo01[which(pela$stage=="Embryo")]<-1 #with embryos
pela<-cbind(pela,Embryo01)
pela$Embryo01<- as.factor(pela$Embryo01)
# 
sapply(datos,class)
#---------------------------------------------------
#Testing different models
#---------------------------------------------------

stag1<-lm(log10(mo2.4)~log10(dry_weight_g),data=pela) # 
summary(stag1)
AIC(stag1)
anova(stag1)
stag2<-lm(log10(mo2.4)~log10(dry_weight_g)*as.numeric(temperature),data=pela)
summary(stag2)
AIC(stag2)
anova(stag2)
stag3<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa.mean,data=pela)
summary(stag3)
AIC(stag3)
anova(stag3)
stag4<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa.mean+as.numeric(temperature),data=pela)
summary(stag4)
anova(stag4)
AIC(stag4)


# Ajusto mejor modelo stag4, ahora considerando huevos como una grupo separado a las otras etapas
eggs<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa.mean+as.numeric(temperature)+Embryo01,data=pela)
summary(eggs)
anova(eggs)
AIC(eggs)
eggs1<-lm(log10(mo2.4)~log10(dry_weight_g)*as.factor(kpa.mean)+as.numeric(temperature)+Embryo01,data=pela)
summary(eggs1)
anova(eggs1)
AIC(eggs1)

AIC(stag4)-AIC(eggs1) # diferencia
# coeficientes y anova
summary(eggs)
Anova(eggs)
# exploratory graphics
plot(log10(mo2.4)~log10(dry_weight_g), data=pela)
# pendiente stag4: todas las etapas como un solo grupo
abline(stag4$coefficients[1]+stag4$coefficients[3]*mean(pela$kpa)+stag4$coefficients[4]*mean(as.numeric(pela$temperature))
       +stag4$coefficients[5]*mean(pela$kpa)*mean(as.numeric(pela$temperature)),stag4$coefficients[2])
# pendiente eggs: huevos como un grupo diferente de las otras tres etapas
abline(eggs$coefficients[1]+eggs$coefficients[3]*mean(pela$kpa)+eggs$coefficients[4]*mean(as.numeric(pela$temperature))
       +eggs$coefficients[6]*mean(pela$kpa)*mean(as.numeric(pela$temperature)),eggs$coefficients[2], col="gray")

#---------------------------------------------------
#plots + predict
setwd("D:/Dropbox/FL+WB/PELA/manuscript")
#setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/Submission")
#setwd("D:/felix.leiva@ulagos.cl/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
{png(filename="1.1. Figure 2.png",width=8,height=4,units="in",res=600)
par(mfrow=c(1,2), tcl=-0.3, family="Arial", oma=c(2,2,0,0))
par(mai=c(0.4,0.5,0.2,0.2))
# PLOT 1

plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
     yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela,!stage=="Embryo"), col="black",cex=0.8)
points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela,stage=="Embryo"), col="gray50",cex=0.8) #display embryo data in grey
# modelo todas las etapas= 1 solo grupo
clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
abline(stag4$coefficients[1]+stag4$coefficients[3]*mean(pela$kpa)+stag4$coefficients[4]*mean(as.numeric(pela$temperature))
       +stag4$coefficients[5]*mean(pela$kpa)*mean(as.numeric(pela$temperature)),stag4$coefficients[2], lwd=2,col="gray50")
# modelo que considera huevos separado de las otras tres etapas
clip(x1=-3.5,x2=0.9,y1=-3,y2=0.95)
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
abline(eggs$coefficients[1]+eggs$coefficients[3]*mean(pela$kpa.mean)+eggs$coefficients[4]*mean(as.numeric(pela$temperature))
       +eggs$coefficients[6]*mean(pela$kpa.mean)*mean(as.numeric(pela$temperature)),eggs$coefficients[2], lwd=2)
# axis labels
mtext(expression(paste(Log[10]," metabolic rate")),side=2,cex=1, outer = T)
mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=T)
# labels slope whit egg
text(-1.5,-3.6,paste("slope =",round(stag4$coefficients[2],2), "?",
                     round(summary(stag4)$coefficients[2,2],2) ," eggs"),cex = 0.8, col="gray50") #Mass exponent=1.05 + 0.016
# # labels slope whithout eggs
text(-1.33,-4,paste("slope =",round(eggs$coefficients[2],2), "?",
                    round(summary(eggs)$coefficients[2,2],2)," no eggs"),cex = 0.8, col="black") #Mass exponent=0.92
text(-4.5,0.8,"A)",font = 2)
#---------------------------------------------------------------------------------------------
# PLOT 2
# Ajusto Modelos por cada nivel de kpa y  a temperatura constante=12°C
#---------------------------------------------------------------------------------------------
#x<-seq(from = min(log10(pela$dry_weight_g)), to = max(log10(pela$dry_weight_g)), length.out = 1000)
# Plot 2
plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),yaxs="i",xaxs="i",col="white",
     xlab="", ylab="",las=1, frame.plot=FALSE,data=pela)
points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela,!stage=="Embryo"), col="black",cex=0.8)
points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela,stage=="Embryo"), col="gray50",cex=0.8) #display embryo data in grey
clip(x1=-3.5,x2=0.9,y1=-3,y2=0.95) # limite para abline
# kpa= 21
df21<-subset(pela,temperature==12 & kpa==21.21)
lmdf21<-lm(log10(mo2.4)~log10(dry_weight_g)+Embryo01,data=df21)
abline(lmdf21$coefficients[1]+lmdf21$coefficients[3]*0,lmdf21$coefficients[2]) # if Embryo=0, without eggs
summary(lmdf21)$coeff
# kpa= 14
df14<-subset(pela,temperature==12 & kpa==14.14)
lmdf14<-lm(log10(mo2.4)~log10(dry_weight_g)+Embryo01,data=df14)
abline(lmdf14$coefficients[1]+lmdf14$coefficients[3]*0,lmdf14$coefficients[2], lty=5) #  # if Embryo=0, without eggs
summary(lmdf14)$coeff
# kpa= 9
df9<-subset(pela,temperature==12 & kpa==9.43)
lmdf9<-lm(log10(mo2.4)~log10(dry_weight_g)+Embryo01,data=df9)#  if Embryo=0, without eggs
abline(lmdf9$coefficients[1]+lmdf9$coefficients[3]*0,lmdf9$coefficients[2], lty=2)
summary(lmdf9)$coeff
# kpa= 5
df5<-subset(pela,temperature==12 & kpa==4.71)
lmdf5<-lm(log10(mo2.4)~log10(dry_weight_g)+Embryo01,data=df5)#without embryos
abline(lmdf5$coefficients[1]+lmdf5$coefficients[3]*0,lmdf5$coefficients[2], lty=4)
summary(lmdf5)$coeff
# kpa= 2
df2<-subset(pela,temperature==12 & kpa==2.36)
lmdf2<-lm(log10(mo2.4)~log10(dry_weight_g)+Embryo01,data=df2)#without embryos
abline(lmdf2$coefficients[1]+lmdf2$coefficients[1]*0,lmdf2$coefficients[2], lty=3)
summary(lmdf2)$coeff

# labels slope models without embryos
clip(x1=-5,x2=0.9,y1=-5,y2=0.95) # limite para abline
text(-1,-3,paste("slope 21 kPa =",round(lmdf21$coefficients[2],2), "?",
                 round(summary(lmdf21)$coefficients[2,2],2)),cex = 0.8, col="black") 
lines(x=c(-3.1,-2.7), y=c(-3,-3))

text(-1,-3.25,paste("slope 14 kPa =",round(lmdf14$coefficients[2],2), "?",
                   round(summary(lmdf14)$coefficients[2,2],2)),cex = 0.8, col="black") 
lines(x=c(-3.1,-2.7), y=c(-3.25,-3.25),lty=5)

text(-1,-3.5,paste("slope"," ","9 kPa =",round(lmdf9$coefficients[2],2), "?",
                   round(summary(lmdf9)$coefficients[2,2],2)),cex = 0.8, col="black") 
lines(x=c(-3.1,-2.7), y=c(-3.5,-3.5),lty=2)

text(-1,-3.75,paste("slope"," " ,"5 kPa =",round(lmdf5$coefficients[2],2), "?",
                   round(summary(lmdf5)$coefficients[2,2],2)),cex = 0.8, col="black") 
lines(x=c(-3.1,-2.7), y=c(-3.75,-3.75),lty=4)

text(-1,-4,paste("slope", " ","2 kPa =",round(lmdf2$coefficients[2],2), "?",
                   round(summary(lmdf2)$coefficients[2,2],2)),cex = 0.8, col="black")  
lines(x=c(-3.1,-2.7), y=c(-4,-4),lty=3)
text(-4.5,0.8,"B)",font = 2)
}
dev.off()


#Figure S2 Supplementary information

temp6.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==6))
temp12.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==12))
temp18.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==18))
{
png(filename="Figure S2.png",width=6,height=6,units="in",res=600)
par(mfrow=c(1,1), tcl=-0.3, family="Arial", oma=c(0,0,0,0))
plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
     yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
points(log10(mo2.4)~log10(dry_weight_g), data=subset(pela, !stage=="Embryo"), col="black",cex=0.8)
points(log10(mo2.4)~log10(dry_weight_g), data=subset(pela, stage=="Embryo"), col="gray50",cex=0.8)
mtext(expression(paste(Log[10]," metabolic rate")),side=2,  line=2, outer=F)
mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=F, line=3)
# Without eggs at 6?C
clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
abline(temp6.s$coefficients[1]+temp6.s$coefficients[3]*mean(pela$kpa)+temp6.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
       temp6.s$coefficients[2], lty=1,col="black")
text(-1.0,-3.5,paste("slope  6°C =",round(temp6.s$coefficients[2],2), "±",
                     round(summary(temp6.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
lines(x=c(-2.3,-2.6), y=c(-3.5,-3.5),lty=1, col="black")

# Without eggs at 12?C
clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
abline(temp12.s$coefficients[1]+temp12.s$coefficients[3]*mean(pela$kpa)+temp12.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
       temp12.s$coefficients[2], lty=2,col="black")
text(-1.0,-3.75,paste("slope 12°C =",round(temp12.s$coefficients[2],2), "±",
                   round(summary(temp12.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
lines(x=c(-2.3,-2.6), y=c(-3.75,-3.75),lty=2, col="black")
# Without eggs at 18?C
clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
abline(temp18.s$coefficients[1]+temp18.s$coefficients[3]*mean(pela$kpa)+temp18.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
       temp18.s$coefficients[2], lty=3,col="black")
text(-1.0,-4,paste("slope 18°C =",round(temp18.s$coefficients[2],2), "±",
                   round(summary(temp18.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
lines(x=c(-2.3,-2.6), y=c(-4,-4),lty=3, col="black")
}
dev.off()






# #Figure S2 Supplementary information
# 
# temp6<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==6))
# temp6.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==6))
# temp12<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==12))
# temp12.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==12))
# temp18<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==18))
# temp18.s<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa+Embryo01,data=subset(pela, temperature==18))
# 
# #png(filename="D:/Dropbox/FL+WB/PELA/manuscript/submission/ 1.1. Figure S2 Scaling of metabolic rate by tempeature.png",width=6,height=6,units="in",res=600)
# plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
#      yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
# points(log10(mo2.4)~log10(dry_weight_g), data=subset(pela, !stage=="Embryo"), col="black",cex=0.8)
# points(log10(mo2.4)~log10(dry_weight_g), data=subset(pela, stage=="Embryo"), col="gray50",cex=0.8)
# mtext(expression(paste(Log[10]," metabolic rate")),side=2,  line=2, outer=F)
# mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=F, line=3)
# 
# # Fit for all life stages at 6?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp6$coefficients[1]+temp6$coefficients[3]*mean(pela$kpa)+temp6$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp6$coefficients[2], lty=1,col="gray50")
# 
# # Fit for all life stages at 12?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp12$coefficients[1]+temp12$coefficients[3]*mean(pela$kpa)+temp12$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp12$coefficients[2], lty=3,col="gray50")
# 
# # Fit for all life stages at 18?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp18$coefficients[1]+temp18$coefficients[3]*mean(pela$kpa)+temp18$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp18$coefficients[2], lty=4,col="gray50")
# 
# text(-0.81,-3,paste("slope 6?C =",round(temp6$coefficients[2],2), "?",
#                     round(summary(temp6)$coefficients[2,2],2)), col="gray50") 
# lines(x=c(-2,-2.3), y=c(-3,-3),lty=1, col="gray50")
# text(-0.79,-3.5,paste("slope 12?C =",round(temp12$coefficients[2],2), "?",
#                       round(summary(temp12)$coefficients[2,2],2)), col="gray50") 
# lines(x=c(-2,-2.3), y=c(-3.5,-3.5),lty=3, col="gray50")
# 
# text(-0.79,-4,paste("slope 18?C =",round(temp18$coefficients[2],2), "?",
#                     round(summary(temp18)$coefficients[2,2],2)), col="gray50")
# lines(x=c(-2,-2.3), y=c(-4,-4),lty=4, col="gray50")
# 
# #dev.off()
# 
# 
# 
# #----------------------------------------------------------------------------------------
# # Modelo para cada temperatura version 2
# #----------------------------------------------------------------------------------------
# temp6<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==6))
# temp12<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==12))
# temp18<-lm(log10(mo2.4)~log10(dry_weight_g)*kpa,data=subset(pela, temperature==18))
# 
# png(filename="D:/Dropbox/FL+WB/PELA/manuscript/submission/ 1.1. Figure S2b Scaling of metabolic rate by tempeature.png",width=8,height=8,units="in",res=600)
# par(mfrow=c(2,2))
# 
# #6?C
# plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
#      yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==6 & !stage=="Embryo"), col="black",cex=0.8)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==6 & stage=="Embryo"), col="gray50",cex=0.8)
# mtext(expression(paste(Log[10]," metabolic rate")),side=2,  line=2, outer=F)
# mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=F, line=3)
# # All life stages at 6?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp6$coefficients[1]+temp6$coefficients[3]*mean(pela$kpa)+temp6$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp6$coefficients[2], lty=1,col="gray50")
# # Without eggs at 6?C
# clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp6.s$coefficients[1]+temp6.s$coefficients[3]*mean(pela$kpa)+temp6.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp6.s$coefficients[2], lty=1,col="black")
# # labels at 6?C
# text(-1.1,-3.5,paste("slope 6?C eggs =",round(temp6$coefficients[2],2), "?",
#                      round(summary(temp6)$coefficients[2,2],2)),cex = 0.9, col="gray50") 
# lines(x=c(-3,-3.3), y=c(-3.5,-3.5),lty=1, col="gray50")
# text(-1.0,-4.0,paste("slope 6?C no eggs =",round(temp6.s$coefficients[2],2), "?",
#                      round(summary(temp6.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
# lines(x=c(-3,-3.3), y=c(-4,-4),lty=1, col="black")
# 
# # 12?C
# plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
#      yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==12 & !stage=="Embryo"), col="black",cex=0.8)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==12 & stage=="Embryo"), col="gray50",cex=0.8)
# mtext(expression(paste(Log[10]," metabolic rates")),side=2,  line=2, outer=F)
# mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=F, line=3)
# 
# # All life stages at 12?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp12$coefficients[1]+temp12$coefficients[3]*mean(pela$kpa)+temp12$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp12$coefficients[2], lty=1,col="gray50")
# # Without eggs at 12?C
# clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp12.s$coefficients[1]+temp12.s$coefficients[3]*mean(pela$kpa)+temp12.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp12.s$coefficients[2], lty=1,col="black")
# # labels at 12C?
# text(-1.1,-3.5,paste("slope 12?C eggs =",round(temp12$coefficients[2],2), "?",
#                      round(summary(temp12)$coefficients[2,2],2)),cex = 0.9, col="gray50") 
# lines(x=c(-3,-3.3), y=c(-3.5,-3.5),lty=1, col="gray50")
# text(-1.0,-4,paste("slope 12?C no eggs =",round(temp12.s$coefficients[2],2), "?",
#                    round(summary(temp12.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
# lines(x=c(-3,-3.3), y=c(-4,-4),lty=1, col="black")
# 
# # Without eggs at 12?C
# plot(log10(mo2.4)~log10(dry_weight_g),cex=0.8,xlim=c(-5,1),ylim=c(-5,1),ylab="", xlab="",
#      yaxs="i",xaxs="i",col="white",las=1,frame.plot=FALSE,data=pela)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==18 & !stage=="Embryo"), col="black",cex=0.8)
# points(log10(mo2.4)~log10(dry_weight_g), data =subset(pela, temperature==18 & stage=="Embryo"), col="gray50",cex=0.8)
# mtext(expression(paste(Log[10]," metabolic rates")),side=2,  line=2, outer=F)
# mtext(expression(paste(Log[10]," body mass")),side=1, cex=1, outer=F, line=3)
# #
# # All life stages at 18?C
# clip(x1=-4.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp18$coefficients[1]+temp18$coefficients[3]*mean(pela$kpa)+temp18$coefficients[4]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp18$coefficients[2], lty=1,col="gray50")
# # Without eggs at 18?C
# clip(x1=-3.5,x2=0.9,y1=-4.5,y2=0.9) # limite para abline
# abline(temp18.s$coefficients[1]+temp18.s$coefficients[3]*mean(pela$kpa)+temp18.s$coefficients[5]*mean(log10(pela$dry_weight_g))*mean(pela$kpa),
#        temp18.s$coefficients[2], lty=1,col="black")
# # labels at 18?C
# text(-1.1,-3.5,paste("slope 18?C eggs =",round(temp18$coefficients[2],2), "?",
#                      round(summary(temp18)$coefficients[2,2],2)),cex = 0.9, col="gray50") 
# lines(x=c(-3,-3.3), y=c(-3.5,-3.5),lty=1, col="gray50")
# text(-1.0,-4,paste("slope 18?C no eggs =",round(temp18.s$coefficients[2],2), "?",
#                    round(summary(temp18.s)$coefficients[2,2],2)),cex = 0.9, col="black") 
# lines(x=c(-3,-3.3), y=c(-4,-4),lty=1, col="black")
# dev.off()
# 

