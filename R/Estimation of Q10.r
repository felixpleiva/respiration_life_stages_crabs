#---------
#Title
#---------
#Respiratory control and thermal sensibility in benthic 
#life stages of porcelain crab (Petrolisthes laevigatus).
#-----------------------------------------------

#---------
#Cleaning working space
rm(list=ls()) 
#Creator: 20160817 FPL
#Modifications: 20170408 CG
today<-format(Sys.Date(),"%Y%m%d")
#--------------------------------------------
#From Linux
# get working directory
setwd("D:/Dropbox/FL+WB/PELA/data")
# Cristóbal
#setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
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
#---------------------------------------------------------------------------

#--------------------------------------------------------
#Transformation of variables wrongly assigned to numeric
#--------------------------------------------------------
pela$stage=as.factor(pela$stage)
pela$oxygen_mg=as.numeric(pela$oxygen_mg)
pela$saturation=as.numeric(pela$saturation)
pela$mo2=as.numeric(pela$mo2)
pela$temperature=as.factor(pela$temperature)
#--------------------------------------------------------

#--------------------------------------------------------
# Mean and standard deviation 
#--------------------------------------------------------
mean1<-aggregate(cbind(mo2,mo2.1, mo2.2,mo2.3,mo2.4)~stage+oxygen_mg+saturation+kpa+temperature,data=pela, mean)
ds1<-aggregate(cbind(mo2,mo2.1, mo2.2,mo2.3,mo2.4)~stage+oxygen_mg+saturation+kpa+temperature,data=pela, sd)
datos<-merge(mean1,ds1,by=c("stage","oxygen_mg","saturation","kpa","temperature"))
names(datos)[6:15]<-c("mo2","mo2.1","mo2.2","mo2.3","mo2.4","mo2.sd","mo2.1.sd",
                      "mo2.2.sd","mo2.3.sd","mo2.4.sd")
#--------------------------------------------------------

#--------------------------------------------------------
# subset by stage
#--------------------------------------------------------
unique(pela$stage)
embrio<-subset(datos,datos$stage=="Embryo")
adult<-subset(datos,datos$stage=="Adult")
juvenil<-subset(datos,datos$stage=="Juvenile")
megalopa<-subset(datos,datos$stage=="Megalopae")
#--------------------------------------------------------

#--------------------------------------------------------
# Calculation of Q10 =(MR2/MR1)^(10/(T2-T1))
#--------------------------------------------------------
# Eggs
#--------------------------------------------------------
embq10<-with(embrio,data.frame(t6.12= (mo2.4[temperature==12]/mo2.4[temperature==6])^(10/(12-6)),
                               t6.18= (mo2.4[temperature==18]/mo2.4[temperature==6])^(10/(18-6)),
                               t12.18=(mo2.4[temperature==18]/mo2.4[temperature==12])^(10/(18-12)),
                               kpa6=kpa[temperature==6],
                               kpa12=kpa[temperature==12],
                               kpa18=kpa[temperature==18])
)
embq10$kpa.mean<- rowMeans(embq10[,c(4,5,6)], na.rm=TRUE)

#--------------------------------------------------------
#Megalopae
#--------------------------------------------------------
megaq10<-with(megalopa,data.frame(t6.12= (mo2.4[temperature==12]/mo2.4[temperature==6])^(10/(12-6)),
                                  t6.18= (mo2.4[temperature==18]/mo2.4[temperature==6])^(10/(18-6)),
                                  t12.18=(mo2.4[temperature==18]/mo2.4[temperature==12])^(10/(18-12)),
                                  kpa6=kpa[temperature==6],
                                  kpa12=kpa[temperature==12],
                                  kpa18=kpa[temperature==18])
)
megaq10$kpa.mean<- rowMeans(megaq10[,c(4,5,6)], na.rm=T)

#--------------------------------------------------------
#Juveniles
#--------------------------------------------------------
juvq10<-with(juvenil,data.frame(t6.12= (mo2.4[temperature==12]/mo2.4[temperature==6])^(10/(12-6)),
                                t6.18= (mo2.4[temperature==18]/mo2.4[temperature==6])^(10/(18-6)),
                                t12.18=(mo2.4[temperature==18]/mo2.4[temperature==12])^(10/(18-12)),
                                kpa6=kpa[temperature==6],
                                kpa12=kpa[temperature==12],
                                kpa18=kpa[temperature==18])
)
juvq10$kpa.mean<- rowMeans(juvq10[,c(4,5,6)], na.rm=TRUE)

#--------------------------------------------------------
# Adults at each oxygen level
#--------------------------------------------------------
aduq10<-with(adult,data.frame(t6.12= (mo2.4[temperature==12]/mo2.4[temperature==6])^(10/(12-6)),
                              t6.18= (mo2.4[temperature==18]/mo2.4[temperature==6])^(10/(18-6)),
                              t12.18=(mo2.4[temperature==18]/mo2.4[temperature==12])^(10/(18-12)), 
                              kpa6=kpa[temperature==6],
                              kpa12=kpa[temperature==12],
                              kpa18=kpa[temperature==18])
)
aduq10$kpa.mean<- rowMeans(aduq10[,c(4,5,6)], na.rm=T)

#---------------------------------------------------------------------------------------------
# Q10 by stages, oxygen (kpa.mean) and delta temperature (DTemp)
#---------------------------------------------------------------------------------------------
Q10<-data.frame(stage= rep(c("Embryo", "Megalopae", "Juvenile", "Adult"),1,each=15),
                DTemp=rep(c(names(embq10)[1:3]),4,each=5),   
                kpa.mean= rep(c(2.36, 4.72, 9.4, 14.2, 21.2),12,each=1),
                Q10= c(embq10[1:5,1], embq10[1:5,2], embq10[1:5,3],
                       megaq10[1:5,1], megaq10[1:5,2] ,megaq10[1:5,3],
                       juvq10[1:5,1], juvq10[1:5,2], juvq10[1:5,3],
                       aduq10[1:5,1], aduq10[1:5,2], aduq10[1:5,3])
)

Q10$stage<-as.factor(Q10$stage)
Q10$DTemp<-as.factor(Q10$DTemp)
#Q10$kpa.mean<-as.factor(Q10$kpa.mean)
#Full model
q10<-subset(Q10, !round(Q10,1)==20.4)
lcu<-with(q10,powerTransform(Q10~kpa.mean*DTemp,family="bcPower")$lambda)
q10=cbind(q10,bcPower(q10$Q10,lcu))
names(q10)[5]<-"bc"
# modelos
###############
lm1.1<-lm(bc~as.factor(kpa.mean)*stage, data=q10)
summary(lm1.1)
Anova(lm1.1)
AIC(lm1.1)
################
lm1<-lm(bc~kpa.mean*stage, data=q10)
summary(lm1)
Anova(lm1)
AIC(lm1)
################
lm.r<-lm(bc~kpa.mean*stage*DTemp, data=q10)
summary(lm.r)
AIC(lm1.1,lm1,lm.r)
Anova(lm.r)
# para respuesta a revisor
res1<-residuals(lm.r, type="pearson")
qqPlot(res1); shapiro.test(res1)
leveneTest(res1~stage*DTemp,data=q10)
# Los datos presentan distribución normal, el unico ruido lo provoca el outlier Q10=20 en adulto
# observo  dentro de etapas
lm.e<-lm(bc~kpa.mean*DTemp, data=subset(q10,stage=="Embryo")); 
lm.m<-lm(bc~kpa.mean*DTemp, data=subset(q10,stage=="Megalopae")); 
lm.j<-lm(bc~kpa.mean*DTemp, data=subset(q10,stage=="Juvenile"));
lm.a<-lm(bc~kpa.mean*DTemp, data=subset(q10,stage=="Adult")); 
# anova embryos
res.e<-residuals(lm.e, type="pearson")
qqPlot(res.e); shapiro.test(res.e)
leveneTest(res.e~DTemp,data=subset(q10,stage=="Embryo"))
Anova(lm.e)
# anova megalopae
res.m<-residuals(lm.m, type="pearson")
qqPlot(res.m); shapiro.test(res.m)
leveneTest(res.m~DTemp,data=subset(q10,stage=="Megalopae"))
Anova(lm.m)
# anova juvenile
res.j<-residuals(lm.j, type="pearson")
qqPlot(res.j); shapiro.test(res.j)
leveneTest(res.j~DTemp,data=subset(q10,stage=="Juvenile"))
Anova(lm.j)
# anova adultos
res.a<-residuals(lm.a, type="pearson")
qqPlot(res.a); shapiro.test(res.a)
leveneTest(res.a~DTemp,data=subset(q10,stage=="Adult"))
Anova(lm.a)
library(lsmeans)
leastsquare = lsmeans(lm.a,pairwise ~ DTemp:kpa.mean)
leastsquare$contrasts
cld(leastsquare)

# observamos medias y medianas (no sensibles a outliers) por variables
# tensión de oxigeno
aggregate(Q10~kpa.mean,data=Q10,mean)
aggregate(Q10~kpa.mean,data=Q10,median)
# etapa
aggregate(Q10~stage,data=Q10,mean)
aggregate(Q10~stage,data=Q10,median)
#--------------------------------------------------------
#Graphics
#--------------------------------------------------------

################################
#END OF SCRIPT
################################
png(filename="Figure 4.png",width=5.8,height=5.5,units="in",res=600)

par(mfrow=c(2,2), tcl=-0.3, family="Arial", oma=c(2,2.2,0,0),frame.plot=FALSE)
omi=c(0.1,0.1,0,0)
# Embryos
par(mai=c(0.2,0.2,0.2,0.2))
# embryos
plot(t6.12~kpa.mean,data=embq10, pch=21,cex=1.4, xlab="",ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",xaxt="n",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=embq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=embq10,pch=21,bg="gray", cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
legend("topleft","Eggs",bty="n")
#Megalopae
plot(t6.12~kpa.mean,data=megaq10, pch=21, cex=1.4, xlab="",ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",xaxt="n",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=megaq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=megaq10, pch=21, bg="gray",cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
legend("topleft","Megalopae",bty="n")
#Juveniles
plot(t6.12~kpa.mean,data=juvq10, pch=21, cex=1.4, ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=juvq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=juvq10, pch=21, bg="gray",cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
legend("topleft","Juveniles",bty="n")
#Adults
# showing a X axis break  
aduq10$t6.12.2<- ifelse(aduq10$t6.12< 8, aduq10$t6.12+10, aduq10$t6.12) #+1
yat<-pretty(aduq10$t6.12.2)
ylab <- ifelse(yat< 19, yat-10, yat) #-17 Y -1
plot(t6.12.2~kpa.mean,data=aduq10, cex=1.4,yaxt="n", ylim=c(10,24),xlim=c(0,25), yaxs="i",xaxs="i",ylab="", xlab="",frame.plot=FALSE)
axis(2, las=1, at=yat, labels=ylab)
axis.break(axis=2,breakpos=19,bgcol="white",breakcol="black", # 17
           style="slash",brw=0.02)
points(t6.18+10~kpa.mean,data=aduq10,cex=1.4,pch=21,bg="black",xaxt="n", ylim=c(20,30))
points(t12.18+10~kpa.mean,data=aduq10,cex=1.4,pch=21,bg="gray",xaxt="n", ylim=c(20,30))
axis(1, labels=F)
mtext("Oxygen tension (kPa)", side=1, outer=T, at=0.5,cex=1,line=1)
mtext(expression(paste(Q)[10]), side=2, outer=T, at=0.5,cex=1,line = 1)
legend("topleft","Adults",bty="n")
dev.off()
#------------------------------------
# plot version 2
#-----------------------------------
setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/Submission")
#png(filename="3.1.Thermal quotient V2 .png",width=7,height=7,units="in",res=600)

# panel 1  (Embriones)
par(mar=c(3,2,2,1), oma=c(2,2,1,2))
par(fig=c(0,5,6.6,9.9)/10)
plot(t6.12~kpa.mean,data=embq10, pch=21,cex=1.4, xlab="",ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",xaxt="n",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=embq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=embq10,pch=21,bg="gray", cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
mtext(expression(paste( italic("(a)"))), side=2, line=1,at=6, las=2)
mtext(expression(paste(Q)[10]), side=2, las=2, line=2)
# panel 2 (Megalopas)
par(fig=c(0,5,3.3,6.6)/10)
par(new=T)
plot(t6.12~kpa.mean,data=megaq10, pch=21, cex=1.4, xlab="",ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",xaxt="n",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=megaq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=megaq10, pch=21, bg="gray",cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
mtext(expression(paste( italic("(b)"))), side=2, line=1,at=6, las=2)
mtext(expression(paste(Q)[10]), side=2, las=2, line=2)
# panel 3 (juveniles)
par(fig=c(0,5,0,3.3)/10)
par(new=T)
plot(t6.12~kpa.mean,data=juvq10, pch=21, cex=1.4, ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,5),yaxs="i",xaxs="i",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=juvq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=juvq10, pch=21, bg="gray",cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
mtext(expression(paste( italic("(c)"))), side=2, line=1,at=6, las=2)
mtext(expression(paste(Q)[10]), side=2, las=2, line=2)
mtext("Oxygen tension (kPa)", side=1, line=2, las=1)
# panel 4 (Adultos)
par(fig=c(5,10,0,9.9)/10)
par(new=T)
plot(t6.12~kpa.mean,data=aduq10, pch=21, cex=1.4, ylab=expression(paste(Q)[10]),
     xlim=c(0,25), ylim=c(0,25),yaxs="i",xaxs="i",las=1,frame.plot=FALSE)
points(t6.18~kpa.mean,data=aduq10, pch=21,bg="black", cex=1.4, ylim=c(0,5), xaxt="n")
points(t12.18~kpa.mean,data=aduq10, pch=21, bg="gray",cex=1.4, ylim=c(0,5), xaxt="n")
axis(1, labels=F)
mtext(expression(paste( italic("(d)"))), side=2, line=1,at=26, las=2)
mtext(expression(paste(Q)[10]), side=2, line=1.5, las=2)
mtext("Oxygen tension (kPa)", side=1, line=2, las=1)
#dev.off()
