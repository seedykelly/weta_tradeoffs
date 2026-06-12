setwd("*******")

##### READ DATA
rm(list=ls())
mydata<-read.csv("******.csv",header=T)  
attach(mydata)
str(mydata)


library(lattice)
library(car)
#########
head.l<-log(head_l)
head.w<-log(head_w)
ear.x<-log(sqrt(ear))
gonad<-log(testes^(1/3))
pro<-log(pronotum)
fem.1<-log(femur.1)
tib.1<-log(tibia.1)
tar.1<-log(tarsus.1)
fem.2<-log(femur.2)
tib.2<-log(tibia.2)
tar.2<-log(tarsus.2)
fem.3<-log(femur.3)
tib.3<-log(tibia.3)
tar.3<-log(tarsus.3)
#eye.x<-sqrt(eye)

weta<-cbind(cat,head.l,pro,head.w,tib.1,fem.1,tib.2,fem.2,tib.3,fem.3,ear.x,gonad,tar.1, tar.2,tar.3)

#you will need to creat a new column called "morph" in which each male is classified as "eight", "nine", or "ten" based head sizes in Kelly & Adams 2010.
#the factor "morph" can then be used in the statistical model

############
weta.data<-weta[,3:12]
model<-prcomp(weta.data,scale=TRUE,cor=TRUE)
model
summary(model)

yv<-predict(model)[,1]
yv2<-predict(model)[,2]
yv3<-predict(model)[,3]
yv4<-predict(model)[,4]
yv5<-predict(model)[,5]
#yv6<-predict(model)[,6]

model.2<-lm(head.l~cat*yv*yv2*yv3*yv4*yv5)
model.2<-lm(head.l~cat+yv+yv2+yv3+yv4+yv5)
summary(model.2)

model.3<-lm(head.l~cat+yv+yv2+yv3+cat:yv:yv3)#
summary(model.3)
anova(model.2,model.3)#

anova(model.3)


plot(yv3~head.l)
abline(lm(yv3~head.l))

############8th instar males
data.2<-cbind(weta,yv)
eight<-subset(data.2,cat=="eight")
weta.8<-eight[,c(2,13)]

model.1<-prcomp(weta.8,scale=TRUE,cor=TRUE)
model.1
summary(model.1)
rel.head8<-predict(model.1)[,2]
x<-data.frame(eight)

ear.1<-lm(ear.x~yv, data=x)
rel.ear8<-ear.1$resid
summary(ear.1)
pro.1<-lm(pro~yv, data=x)
rel.pro8<-pro.1$resid

cor(rel.ear8,head.l)
plot(rel.ear8~head.l,zzz)
abline(lm(rel.ear8~head.l,zzz))
ccc<-lm(rel.ear8~head.l,zzz)
summary(ccc)


###########9th instar males
data.2<-cbind(weta,yv)
nine<-subset(data.2,cat=="nine")
weta.9<-nine[,c(2,13)]
cor(nine)
model.1<-prcomp(weta.9,scale=TRUE,cor=TRUE)
model.1
summary(model.1)
rel.head9<-predict(model.1)[,2]
y<-data.frame(nine)

ear.1<-lm(ear.x~yv, data=y)
rel.ear9<-ear.1$resid
summary(ear.1)
pro.1<-lm(pro~yv, data=y)
rel.pro9<-pro.1$resid


######10th instar
data.2<-cbind(weta,yv)
ten<-subset(data.2,cat=="ten")
weta.10<-ten[,c(2,13)]
cor(ten)
model.1<-prcomp(weta.10,scale=TRUE,cor=TRUE)
model.1
summary(model.1)
rel.head10<-predict(model.1)[,2]
z<-data.frame(ten)

ear.1<-lm(ear.x~yv, data=z)
rel.ear10<-ear.1$resid

pro.1<-lm(pro~yv, data=z)
rel.pro10<-pro.1$resid


