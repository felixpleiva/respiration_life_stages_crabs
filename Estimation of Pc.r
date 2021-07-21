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
# setwd("~/Insync/Laboratorio.Poblaciones.Marinas/Publicaciones/Plasticidad pela/analisis/datos")
getwd()     #verifico directorio
#----------------------------------------------

#-----------------------------------------------
#Libraries
#-----------------------------------------------
library(plotrix);library(ggplot2);library(nlme);
library(car);library(multcomp); 
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

pela$kpa.mean<-ifelse(pela$kpa==2.07 | pela$kpa==2.36 | pela$kpa==2.65, 2.36, 
               ifelse(pela$kpa==4.14 | pela$kpa==4.71 | pela$kpa==5.30, 4.71,
               ifelse(pela$kpa==8.29 | pela$kpa==9.43 | pela$kpa==10.60, 9.44,
               ifelse(pela$kpa==12.43 | pela$kpa==14.14 | pela$kpa==15.89, 14.15,
               ifelse(pela$kpa==18.65 | pela$kpa==21.21 | pela$kpa==23.84, 21.23, NA)))))

table(pela$stage, pela$kpa.mean)
# ordeno niveles etapa
pela$stage<- factor(pela$stage, levels = c("Embryo","Megalopae","Juvenile","Adult"))
#-----------------------------------------------------------------------
# Estimo Pc y creo graficos 
#-----------------------------------------------------------------------------
library(segmented)
# Adultos
# Temperatura = 6
t6<-subset(pela, temperature==6 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t6)
# modelo lineal
lm6.1<-lm(mo2.3~kpa.mean,data=t6)
abline(lm6.1$coefficients,col="red")
# polinomial
lm6.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t6)
lines(sort(t6$kpa.mean), fitted(lm6.2)[order(t6$kpa.mean)], col='blue', type='b') 
#
AIC(lm6.1, lm6.2)
# Estima P critico
pcrit1 <- segmented(lm6.1, 
                    seg.Z = ~ kpa.mean, 
                    psi = 9)
summary(pcrit1)
pcrit1$psi
pcrit2 <- segmented(lm6.2, 
                    seg.Z = ~ poly(kpa.mean, 3), 
                    psi = 4)
abline(v=pcrit1$psi[2], lwd=2) # Pcritico???


# Adultos
# Temperatura = 12
t12<-subset(pela, temperature==12 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t12)
# modelo lineal
lm12.1<-lm(mo2.3~kpa.mean,data=t12)
abline(lm12.1$coefficients,col="red")
# polinomial
lm12.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t12)
lines(sort(t12$kpa.mean), fitted(lm12.2)[order(t12$kpa.mean)], col='blue', type='b') 
#
AIC(lm12.1, lm12.2)
# Estima P critico
pcrit_t12a <- segmented(lm12.1, 
                    seg.Z = ~ kpa.mean, 
                    psi =7)
summary(pcrit_t12a)
pcrit_t12a$psi


# Adultos
# Temperatura = 18
t18<-subset(pela, temperature==18 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t18)
# modelo lineal
lm18.1<-lm(mo2.3~kpa.mean,data=t18)
abline(lm18.1$coefficients,col="red")
# polinomial
lm18.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t18)
lines(sort(t18$kpa.mean), fitted(lm18.2)[order(t18$kpa.mean)], col='blue', type='b') 
#
AIC(lm18.1, lm18.2)
# Estima P critico
pcrit_t18a <- segmented(lm18.1, 
                        seg.Z = ~ kpa.mean, 
                        psi =14)
summary(pcrit_t12a)
pcrit_t12a$psi   
#---------------------------------------------------
## Trabajo con media por tratamiento
#--------------------------------------------
pela1<-aggregate(cbind(mo2,mo2.1, mo2.2,mo2.3,mo2.4)~stage+kpa.mean+temperature,data=pela, mean)
# Adultos
# Temperatura = 6
t6<-subset(pela1, temperature==6 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t6)
# modelo lineal
lm6.1<-lm(mo2.3~kpa.mean,data=t6)
abline(lm6.1$coefficients,col="red")
# polinomial
lm6.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t6)
lines(sort(t6$kpa.mean), fitted(lm6.2)[order(t6$kpa.mean)], col='blue', type='b') 
#
AIC(lm6.1, lm6.2)
# Estima P critico
pcrit1 <- segmented(lm6.1, 
                    seg.Z = ~ kpa.mean, 
                    psi = 9)
summary(pcrit1)
pcrit1$psi
# no  resulto con el modelo polinomial
pcrit2 <- segmented(lm6.2, 
                    seg.Z = ~ poly(kpa.mean, 3), 
                    psi = 4)
abline(v=pcrit1$psi[2], lwd=2) # Pcritico???


# Adultos
# Temperatura = 12
t12<-subset(pela1, temperature==12 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t12)
# modelo lineal
lm12.1<-lm(mo2.3~kpa.mean,data=t12)
abline(lm12.1$coefficients,col="red")
# polinomial
lm12.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t12)
lines(sort(t12$kpa.mean), fitted(lm12.2)[order(t12$kpa.mean)], col='blue', type='b') 
#
AIC(lm12.1, lm12.2)
# Estima P critico
pcrit_t12a <- segmented(lm12.1, 
                        seg.Z = ~ kpa.mean, 
                        psi =5)
# Errorrrrrrrrrrrrrrr


# Adultos
# Temperatura = 18
t18<-subset(pela1, temperature==18 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t18)
# modelo lineal
lm18.1<-lm(mo2.3~kpa.mean,data=t18)
abline(lm18.1$coefficients,col="red")
# polinomial
lm18.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t18)
lines(sort(t18$kpa.mean), fitted(lm18.2)[order(t18$kpa.mean)], col='blue', type='b') 
#
AIC(lm18.1, lm18.2)
# Estima P critico
pcrit_t18a <- segmented(lm18.1, 
                        seg.Z = ~ kpa.mean, 
                        psi =8)
# ERROR
# psi  es un valor inicial que se le da

#--------------------------------------------------------







#--------------------------------------------------------------------
## Calculamos Pc de acuerdo a Mueller and Seymor (2011).
#The regulation index: a new method for assessing 
#the relationship between oxygen consumption and environmental oxygen
#--------------------------------------------------------------------
# Adultos
# Temperatura = 6
t6<-subset(pela, temperature==6 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t6)
# modelo lineal
lm6.1<-lm(mo2.3~kpa.mean,data=t6)
#abline(lm6.1$coefficients,col="red")
# linea recta ("Perfect conformity")
v1<-data.frame(aggregate(mo2.3~kpa.mean, t6,mean))[c(1,5),]
m1<-lm(mo2.3~kpa.mean,data=v1)
abline(m1, col="blue")
# polinomial
lm6.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t6)
AIC(lm6.1, lm6.2)
lines(sort(t6$kpa.mean), fitted(lm6.2)[order(t6$kpa.mean)], col='blue', type='b') 
# Estima P critico (Distancia maxima entre linea recta y la curva)
pre.t6=data.frame(kpa.mean=seq(2.36,21.23,0.1))
pre.t6$curva<-predict(lm6.2,newdata = pre.t6)
pre.t6$recta<-predict(m1,newdata = pre.t6)
lines(pre.t6$kpa.mean,pre.t6$curva, col="red")
lines(pre.t6$kpa.mean,pre.t6$recta, col="red")
# observo el punto de mayor diferencia 
pre.t6$dif<-pre.t6$curva-pre.t6$recta
subset(pre.t6,dif==max(pre.t6$dif)) # P.critico= 9.16 kpa

# Adultos
# Temperatura = 12
t12<-subset(pela, temperature==12 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t12)
# modelo lineal
lm12.1<-lm(mo2.3~kpa.mean,data=t12)
#abline(lm6.1$coefficients,col="red")
# linea recta ("Perfect conformity")
v2<-data.frame(aggregate(mo2.3~kpa.mean, t12,mean))[c(1,5),]
m2<-lm(mo2.3~kpa.mean,data=v2)
abline(m2, col="blue")
# polinomial
lm12.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t12)
AIC(lm12.1, lm12.2)
lines(sort(t12$kpa.mean), fitted(lm12.2)[order(t12$kpa.mean)], col='blue', type='b') 
# Estima P critico (Distancia maxima entre linea recta y la curva)
pre.t12=data.frame(kpa.mean=seq(2.36,21.23,0.1))
pre.t12$curva<-predict(lm12.2,newdata = pre.t12)
pre.t12$recta<-predict(m2,newdata = pre.t12)
lines(pre.t12$kpa.mean,pre.t12$curva, col="red")
lines(pre.t12$kpa.mean,pre.t12$recta, col="red")
# observo el punto de mayor diferencia 
pre.t12$dif<-pre.t12$curva-pre.t12$recta
subset(pre.t12,dif==max(pre.t12$dif)) # P.critico= 15.46 kpa

# Adultos
# Temperatura = 18
t18<-subset(pela, temperature==18 & stage=="Adult")
plot(mo2.3~kpa.mean,data=t18)
# modelo lineal
lm18.1<-lm(mo2.3~kpa.mean,data=t18)
#abline(lm6.1$coefficients,col="red")
# linea recta ("Perfect conformity")
v3<-data.frame(aggregate(mo2.3~kpa.mean, t18,mean))[c(1,5),]
m3<-lm(mo2.3~kpa.mean,data=v3)
abline(m3, col="blue")
# polinomial
lm18.2<-lm(mo2.3 ~ poly(kpa.mean, 3), data = t18)
AIC(lm18.1, lm18.2)
lines(sort(t18$kpa.mean), fitted(lm18.2)[order(t18$kpa.mean)], col='blue', type='b') 
# Estima P critico (Distancia maxima entre linea recta y la curva)
pre.t18=data.frame(kpa.mean=seq(2.36,21.23,0.1))
pre.t18$curva<-predict(lm18.2,newdata = pre.t18)
pre.t18$recta<-predict(m3,newdata = pre.t18)
lines(pre.t18$kpa.mean,pre.t18$curva, col="red")
lines(pre.t18$kpa.mean,pre.t18$recta, col="red")
# observo el punto de mayor diferencia 
pre.t18$dif<-pre.t18$curva-pre.t18$recta
subset(pre.t18,dif==max(pre.t18$dif)) # P.critico= 15.76 kpa
