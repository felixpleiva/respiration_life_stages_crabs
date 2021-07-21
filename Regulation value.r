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
#Modifications: 20170423 CG
today<-format(Sys.Date(),"%Y%m%d")
#--------------------------------------------
#From Linux
# get working directory
#setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
# getwd()     #verifico directorio
#----------------------------------------------
setwd("D:/Dropbox/FL+WB/PELA/data")
getwd()     #verifico directorio

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


# porcentaje mo2.3 por temperatura [embriones]
embrio$mo2_max<- ifelse(embrio$temperature==6,100*embrio$mo2.3/max(embrio$mo2.3[embrio$temperature==6]),
                 ifelse(embrio$temperature==12,100*embrio$mo2.3/max(embrio$mo2.3[embrio$temperature==12]),
                 ifelse(embrio$temperature==18,100*embrio$mo2.3/max(embrio$mo2.3[embrio$temperature==18]),NA)))
# porcentaje mo2.3 por temperatura [megalopa]                 
megalopa$mo2_max<- ifelse(megalopa$temperature==6,100*megalopa$mo2.3/max(megalopa$mo2.3[megalopa$temperature==6]),
                   ifelse(megalopa$temperature==12,100*megalopa$mo2.3/max(megalopa$mo2.3[megalopa$temperature==12]),
                   ifelse(megalopa$temperature==18,100*megalopa$mo2.3/max(megalopa$mo2.3[megalopa$temperature==18]),NA)))
# porcentaje mo2.3 por temperatura [juvenile]                 
juvenil$mo2_max<- ifelse(juvenil$temperature==6,100*juvenil$mo2.3/max(juvenil$mo2.3[juvenil$temperature==6]),
                  ifelse(juvenil$temperature==12,100*juvenil$mo2.3/max(juvenil$mo2.3[juvenil$temperature==12]),
                  ifelse(juvenil$temperature==18,100*juvenil$mo2.3/max(juvenil$mo2.3[juvenil$temperature==18]),NA)))
# porcentaje mo2.3 por temperatura [adulto]  
adult$mo2_max<- ifelse(adult$temperature==6,100*adult$mo2.3/max(adult$mo2.3[adult$temperature==6]),
                ifelse(adult$temperature==12,100*adult$mo2.3/max(adult$mo2.3[adult$temperature==12]),
                ifelse(adult$temperature==18,100*adult$mo2.3/max(adult$mo2.3[adult$temperature==18]),NA)))
#############################################################################################
# Fija modelo no lineal y calcula area bajo la curva
############################################################################################
#jpeg(paste("2.1 Figure 2 Oxygen regulation value.jpeg"),height = 18,width = 18,units="cm",quality=100,res=300)
par(mfrow=c(2,2),oma=c(2,3,1,1),mar=c(2,1,0,0)) #mar fija margenes alrededor de figuras
# # #1.0.Embryos
#----------------------------------------------------------------------------------------
#Temperature= 6
#m2_max=A*saturation^3 + B*sauration^2 + C*saturation
x6=data.frame(saturation=seq(0,max(embrio$saturation[embrio$temperature==6]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(embrio,temperature==6),start=list(A=10,B=1,C=0))
summary(mm2)
y6=c(predict(mm2,newdata=x6))
#area bajo la curva
R.em.6<-sum(y6)/(dim(x6)[1]*100)*100
#Temperature= 12
x12=data.frame(saturation=seq(0,max(embrio$saturation[embrio$temperature==12]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(embrio,temperature==12),start=list(A=10,B=1,C=0))
summary(mm2)
y12=c(predict(mm2,newdata=x12))
#area bajo la curva
R.em.12<-sum(y12)/(dim(x12)[1]*100)*100
#Temperature= 18
x18=data.frame(saturation=seq(0,max(embrio$saturation[embrio$temperature==18]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(embrio,temperature==18),start=list(A=10,B=1,C=0))
summary(mm2)
y18=c(predict(mm2,newdata=x18))
#area bajo la curva
R.em.18<-sum(y18)/(dim(x18)[1]*100)*100  
#---------------------------------------------------------------------------------------------
#PLOT
#----------------------------------------------------------------------------------------------
# defino color gris transparente (cambia color del alpha)
col2rgb("gray", alpha=TRUE) # 
graytrans <- rgb(190, 190, 190, 127, maxColorValue=255) #doy transparencia al relleno del polygono
# plot [t=6]
with(subset(embrio,embrio$temperature==6), 
      plot(saturation,mo2_max,pch=21,cex=1.2,bg="blue",xaxs="i", yaxs="i",
      ylim=c(0,102),xlim=c(0,102),las=1,
      ylab="",
      xlab="",xaxt="n",
       polygon(c(x6[,1],max(x6)),c(y6,0),col=graytrans, border="blue",lty=2))
     )

# plot [t=12]
with(subset(embrio,embrio$temperature==12), 
      points(saturation,mo2_max,pch=21,bg="yellow",cex=1.2,xaxs="i", yaxs="i",
             ylim=c(0,102),xlim=c(0,102),las=1,
             ylab="",
             polygon(c(x12[,1],max(x12)),c(y12,0),col=graytrans,border="yellow",lty=2)))
 
#plot [t=18]
with(subset(embrio,embrio$temperature==18), 
      points(saturation,mo2_max,pch=21,bg="red",cex=1.2,xaxs="i", yaxs="i",
             ylim=c(0,102),xlim=c(0,102),las=1,
             ylab="",
             polygon(c(x18[,1],max(x18)),c(y18,0),col=graytrans,border="red",lty=2))
     )
text(10,95,paste("R=",round(R.em.6,2),"%",sep=""),cex = 1, col="blue")
text(10,90,paste("R=",round(R.em.12,2),"%",sep=""),cex = 1, col="yellow")
text(10,85,paste("R=",round(R.em.18,2),"%",sep=""),cex = 1, col="red")
# dev.off()
# #--------------------------------------
# #2.0.Megalopae
# #--------------------------------------
#Temperature= 6
#m2_max=A*saturation^3 + B*sauration^2 + C*saturation
x6=data.frame(saturation=seq(0,max(megalopa$saturation[megalopa$temperature==6]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(megalopa,temperature==6),start=list(A=10,B=1,C=0))
summary(mm2)
y6=c(predict(mm2,newdata=x6))
#area bajo la curva
R.meg.6<-sum(y6)/(dim(x6)[1]*100)*100
#Temperature= 12
x12=data.frame(saturation=seq(0,max(megalopa$saturation[megalopa$temperature==12]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(megalopa,temperature==12),start=list(A=10,B=1,C=0))
summary(mm2)
y12=c(predict(mm2,newdata=x12))
#area bajo la curva
R.meg.12<-sum(y12)/(dim(x12)[1]*100)*100
#Temperature= 18
x18=data.frame(saturation=seq(0,max(megalopa$saturation[megalopa$temperature==18]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(megalopa,temperature==18),start=list(A=10,B=1,C=0))
summary(mm2)
y18=c(predict(mm2,newdata=x18))
#area bajo la curva
R.meg.18<-sum(y18)/(dim(x18)[1]*100)*100  
#-------------------------------------
#PLOT
#--------------------------------------
# plot [t=6]
with(subset(megalopa,megalopa$temperature==6), 
     plot(saturation,mo2_max,pch=21,cex=1.2,bg="blue",xaxs="i", yaxs="i",
          ylim=c(0,102),xlim=c(0,102),las=1,
          ylab="", yaxt="n",
          xlab="", xaxt="n",
          polygon(c(x6[,1],max(x6)),c(y6,0),col=graytrans, border="blue",lty=2))
)

# plot [t=12]
with(subset(megalopa,megalopa$temperature==12), 
     points(saturation,mo2_max,pch=21,bg="yellow",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x12[,1],max(x12)),c(y12,0),col=graytrans,border="yellow",lty=2)))

#plot [t=18]
with(subset(megalopa,megalopa$temperature==18), 
     points(saturation,mo2_max,pch=21,bg="red",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x18[,1],max(x18)),c(y18,0),col=graytrans,border="red",lty=2))
)
text(10,95,paste("R=",round(R.meg.6,2),"%",sep=""),cex = 1, col="blue")
text(9,90,paste("R=",round(R.meg.12,2),"%",sep=""),cex = 1, col="yellow")
text(10,85,paste("R=",round(R.meg.18,2),"%",sep=""),cex = 1, col="red")
# #--------------------------------------
# #3.0.Juveniles
# #--------------------------------------
#Temperature= 6
#m2_max=A*saturation^3 + B*sauration^2 + C*saturation
x6=data.frame(saturation=seq(0,max(juvenil$saturation[juvenil$temperature==6]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(juvenil,temperature==6),start=list(A=10,B=1,C=0))
summary(mm2)
y6=c(predict(mm2,newdata=x6))
#area bajo la curva
R.juv.6<-sum(y6)/(dim(x6)[1]*100)*100
#Temperature= 12
x12=data.frame(saturation=seq(0,max(juvenil$saturation[juvenil$temperature==12]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(juvenil,temperature==12),start=list(A=10,B=1,C=0))
summary(mm2)
y12=c(predict(mm2,newdata=x12))
#area bajo la curva
R.juv.12<-sum(y12)/(dim(x12)[1]*100)*100
#Temperature= 18
x18=data.frame(saturation=seq(0,max(juvenil$saturation[juvenil$temperature==18]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(juvenil,temperature==18),start=list(A=10,B=1,C=0))
summary(mm2)
y18=c(predict(mm2,newdata=x18))
#area bajo la curva
R.juv.18<-sum(y18)/(dim(x18)[1]*100)*100  
#-------------------------------------
#PLOT
#--------------------------------------
# plot [t=6]
with(subset(juvenil,juvenil$temperature==6), 
     plot(saturation,mo2_max,pch=21,cex=1.2,bg="blue",xaxs="i", yaxs="i",
          ylim=c(0,102),xlim=c(0,102),las=1,
          ylab="",
          xlab="", 
          polygon(c(x6[,1],max(x6)),c(y6,0),col=graytrans, border="blue",lty=2))
)

# plot [t=12]
with(subset(juvenil,juvenil$temperature==12), 
     points(saturation,mo2_max,pch=21,bg="yellow",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x12[,1],max(x12)),c(y12,0),col=graytrans,border="yellow",lty=2)))

#plot [t=18]
with(subset(juvenil,juvenil$temperature==18), 
     points(saturation,mo2_max,pch=21,bg="red",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x18[,1],max(x18)),c(y18,0),col=graytrans,border="red",lty=2))
)
text(10,95,paste("R=",round(R.juv.6,2),"%",sep=""),cex = 1, col="blue")
text(10,90,paste("R=",round(R.juv.12,2),"%",sep=""),cex = 1, col="yellow")
text(10,85,paste("R=",round(R.juv.18,2),"%",sep=""),cex = 1, col="red")
# #--------------------------------------
# #4.0.Adults
# #--------------------------------------
#Temperature= 6
#m2_max=A*saturation^3 + B*sauration^2 + C*saturation
x6=data.frame(saturation=seq(0,max(adult$saturation[adult$temperature==6]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(adult,temperature==6),start=list(A=10,B=1,C=0))
summary(mm2)
y6=c(predict(mm2,newdata=x6))
#area bajo la curva
R.adu.6<-sum(y6)/(dim(x6)[1]*100)*100
#Temperature= 12
x12=data.frame(saturation=seq(0,max(adult$saturation[adult$temperature==12]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(adult,temperature==12),start=list(A=10,B=1,C=0))
summary(mm2)
y12=c(predict(mm2,newdata=x12))
#area bajo la curva
R.adu.12<-sum(y12)/(dim(x12)[1]*100)*100
#Temperature= 18
x18=data.frame(saturation=seq(0,max(adult$saturation[adult$temperature==18]),0.001))
mm2=nls(mo2_max~A*saturation^3+B*saturation^2+C*saturation,
        data=subset(adult,temperature==18),start=list(A=10,B=1,C=0))
summary(mm2)
y18=c(predict(mm2,newdata=x18))
#area bajo la curva
R.adu.18<-sum(y18)/(dim(x18)[1]*100)*100  
#-------------------------------------
#PLOT
#-------------------------------------
# plot [t=6]
with(subset(adult,adult$temperature==6), 
     plot(saturation,mo2_max,pch=21,cex=1.2,bg="blue",xaxs="i", yaxs="i",
          ylim=c(0,102),xlim=c(0,102),las=1,
          ylab="", yaxt="n",
          xlab="",
          polygon(c(x6[,1],max(x6)),c(y6,0),col=graytrans, border="blue",lty=2))
)

# plot [t=12]
with(subset(adult,adult$temperature==12), 
     points(saturation,mo2_max,pch=21,bg="yellow",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x12[,1],max(x12)),c(y12,0),col=graytrans,border="yellow",lty=2)))

#plot [t=18]
with(subset(adult,adult$temperature==18), 
     points(saturation,mo2_max,pch=21,bg="red",cex=1.2,xaxs="i", yaxs="i",
            ylim=c(0,102),xlim=c(0,102),las=1,
            ylab="",
            polygon(c(x18[,1],max(x18)),c(y18,0),col=graytrans,border="red",lty=2))
)
text(10,95,paste("R=",round(R.adu.6,2),"%",sep=""),cex = 1, col="blue")
text(10,90,paste("R=",round(R.adu.12,2),"%",sep=""),cex = 1, col="yellow")
text(10,85,paste("R=",round(R.adu.18,2),"%",sep=""),cex = 1, col="red")
mtext(expression(paste("Maximun ",dot(M)[O[2]],"(%)")),side = 2,line =8.5,adj = 0)
mtext("Oxygen saturation (%)",side = 1,line =2,adj = 1)
#dev.off()

#---------------------------------------------------------------------------------------------
# Regulation value by stage and temperature
#---------------------------------------------------------------------------------------------
RV<-data.frame(stage= rep(c("Embryo", "Megalopae", "Juvenile", "Adult"),1,each=3),
               temperature= rep(c(6,12,18),4,each=1),
               RV= c(R.em.6,R.em.12, R.em.18,
                     R.meg.6, R.meg.12 ,R.meg.18,
                     R.juv.6, R.juv.12, R.juv.18,
                     R.adu.6, R.adu.12, R.adu.18))

lm.R1<-lm(RV~as.factor(temperature)+stage, data=RV) 
summary(lm.R1)
Anova(lm.R1)
lm.R2<-lm(RV~as.numeric(temperature)+stage, data=RV) 
summary(lm.R2)
Anova(lm.R2)
AIC(lm.R1,lm.R2)


res.R<-residuals(lm.R, type="pearson")
library(car)
qqPlot(res.R); shapiro.test(res.R)
leveneTest(res.R~stage,data=RV)
plot(res.R~RV$stage)
#--------------------------------------------------------
# Predicciones
#--------------------------------------------------------
RV$pred<-predict(lm.R,newdata=RV,re.form=NA)
RV<-cbind(RV,data.frame(predict(lm.R, newdata=RV,se.fit=T)[2]))
#--------------------------------------------------------
#Graphic
# #--------------------------------------------------------
setwd("D:/Dropbox/FL+WB/PELA/manuscript")
png(filename="Figure 3v2.png",width=5.8,height=5.8,units="in",res=600)
par(mfrow=c(2,2), tcl=-0.3, family="Arial", oma=c(2,2.2,0,0))
omi=c(0.1,0.1,0,0)
# Embryos
par(mai=c(0.2,0.2,0.2,0.2))
plot(RV~as.numeric(as.character((temperature))),data=subset(RV,stage=="Embryo"), ylim=c(0,100),
     xlim=c(0,20),ylab="", yaxs="i",xaxs="i",xaxt="n",xlab="", las=1,frame.plot=FALSE)
axis(1,at=c(0,5,10,15,20),labels=c("","","","",""))
abline(h=50,lty=5)
text(2,90,"Eggs")
text(15,20,expression(italic("Hypoxia sensitive")), cex=0.9, col="gray30")
text(15,80,expression(italic("Oxygen regulator")), cex=0.9, col="gray30")
# megalopa
plot(RV~as.numeric(as.character((temperature))),data=subset(RV,stage=="Megalopae"), ylim=c(0,100),
     xlim=c(0,20),ylab="",yaxs="i",xaxs="i",xaxt="n",xlab="", las=1,frame.plot=FALSE)
axis(1,at=c(0,5,10,15,20),labels=c("","","","",""))
abline(h=50,lty=5)
text(3.5,90,"Megalopae")
text(15,20,expression(italic("Hypoxia sensitive")), cex=0.9, col="gray30")
text(15,80,expression(italic("Oxygen regulator")), cex=0.9, col="gray30")

# juvenile
plot(RV~as.numeric(as.character((temperature))),data=subset(RV,stage=="Juvenile"), ylim=c(0,100),
     xlim=c(0,20),ylab="",yaxs="i",xaxs="i",xaxt="n",xlab="", las=1,frame.plot=FALSE)
axis(1,at=c(0,5,10,15,20),labels=c(0,5,10,15,20))
abline(h=50,lty=5)
text(3,90,"Juveniles")
text(15,20,expression(italic("Hypoxia sensitive")), cex=0.9, col="gray30")
text(15,80,expression(italic("Oxygen regulator")), cex=0.9, col="gray30")

# adult
plot(RV~as.numeric(as.character((temperature))),data=subset(RV,stage=="Adult"), ylim=c(0,100),
     xlim=c(0,20),ylab="",yaxs="i",xaxs="i",xaxt="n",xlab="", las=1,frame.plot=FALSE)
abline(h=50,lty=5)
text(2,90,"Adults")
axis(1,at=c(0,5,10,15,20),labels=c(0,5,10,15,20))
text(15,20,expression(italic("Hypoxia sensitive")), cex=0.9, col="gray30")
text(15,80,expression(italic("Oxygen regulator")), cex=0.9, col="gray30")

mtext("Temperature (°C)",side=1, outer=T, line=1 ,adj=0.5)
mtext("Regulation value (%)", side=2, outer=T, at=0.5,cex=1,line = 1)
dev.off()



