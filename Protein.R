source('summarySE.R')

#install.packages("MCMCglmm")
library(MCMCglmm)
library(plyr)
#install.packages("car") #takes a while
library(car)
#install.packages("multcomp")
library(multcomp)
#install.packages("gtools")
library(gtools)
#install.packages("quantreg")
library(quantreg)
#install.packages("calibrate")
library(calibrate)
#install.packages("MASS")
library(MASS)
#install.packages("AICcmodavg") #takes a while
library(AICcmodavg)
#install.packages("e1071")
library(e1071)
#install.packages("nlme")
library(nlme)
#install.packages("lmmfit") #incompatible with current version
library(lmmfit)
#install.packages("labdsv")
library(labdsv)
#install.packages("vegan")
library(vegan)
#install.packages("plotrix")
library(plotrix)
#install.packages("pgirmess")
library(pgirmess)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("lme4")
library(lme4)
#install.packages("knitr")
library(knitr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("egg")
library(egg)
#install.packages("MuMIn")
library(MuMIn)

##############AdultProtein##############

AdultPro=read.csv("~/Dropbox/PoritesSpawnMay2019/Protein/AdtPro.csv")
str(AdultPro)
AdultPro$Fam=as.factor(AdultPro$Fam)
AdultPro$Rep=as.factor(AdultPro$Rep)
hist(AdultPro$PerArea)
AdultPro[AdultPro$PerArea>160,] #H3-20 had very high conc
#AdultSym$Log=log(AdultSym$Conc)
#hist(AdultSym$Log) #still doesn't look normally distributed
#AdultSym=AdultSym[-c(56), ]
#hist(AdultSym$Log)

quartz()

g7=ggplot(AdultPro,aes(factor(Trmt),PerArea,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(µg~protein~per~cm^2))+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")+ #remove legend
  theme(axis.title.x = element_blank()) #remove x axis title

all.marg=summarySE(AdultPro,measurevar="PerArea",groupvars=c("Trmt","Origin"), na.rm=T)
all.marg #inshore=(116.83941-118.07787)/118.07787=-0.0104885 offshore=(98.26623-129.13152)/129.13152=-0.2390221 dif=-0.0104885-(-0.2390221)=0.2285336

all.marg1=summarySE(AdultPro,measurevar="PerArea",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg1

pd <- position_dodge(.3)
quartz()

p2=ggplot(all.marg1,aes(x=Fam,y=PerArea,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=PerArea+se,ymin=PerArea-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(µg~protein~per~cm^2),title="Adult")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())

###Stats###
model1<- lme(PerArea~Trmt*Origin,data=AdultPro,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model1) # Sig interaction

# Fixed effects:  PerArea ~ Trmt * Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)             118.07787 10.428755 48 11.322336  0.0000
# TrmtHeat                 -1.23846  9.791917 48 -0.126478  0.8999
# OriginOffshore           11.05365 14.748487  8  0.749477  0.4750
# TrmtHeat:OriginOffshore -30.12624 13.843437 48 -2.176211  0.0345


#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model1,which=1) 

r.squaredGLMM(model1) #Fam: 0.1140402 0.381456 # Fam + Tank: 0.1145011 0.9537757

##############LarvaeProtein##############

LarvaePro=read.csv("~/Dropbox/PoritesSpawnMay2019/Protein/LarPro.csv")
str(LarvaePro)
LarvaePro$Fam=as.factor(LarvaePro$Fam)
LarvaePro$Rep=as.factor(LarvaePro$Rep)

### Pre-exp larvae
LarDRPro=subset(LarvaePro,Trmt==c("DR"),select=c(Fam,Trmt,Origin,Rep,PrPerLar))
hist(LarDRPro$PrPerLar)
LarDRPro$Log=log(LarDRPro$PrPerLar)
hist(LarDRPro$Log)

quartz()

ggplot(LarDRPro,aes(factor(Origin),Log))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Origin",y="Log (µg protein per larvae)")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))

all.marg=summarySE(LarDRPro,measurevar="Log",groupvars=c("Fam","Origin"), na.rm=T)
all.marg

pd <- position_dodge(.2)
quartz()

ggplot(all.marg,aes(x=Fam,y=Log,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(x="Family",y="Log (µg protein per larvae)")+
  theme_bw()+ 
  scale_shape_manual(values = c(2,8))

###Stats###
model<- lme(Log~Origin,data=LarDRPro,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model) 

# Fixed effects:  Log ~ Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)     1.2023166 0.1034293 20 11.624532  0.0000
# OriginOffshore -0.0866664 0.1462710  8 -0.592506  0.5699

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model) #0.02198739 0.4411293


###Experiment
LarExp=subset(LarvaePro,Trmt==c("Control","Heat"),select=c(Fam,Trmt,Origin,Rep,PrPerLar))
hist(LarExp$PrPerLar)
LarExp[LarExp$PrPerLar>6,] #Fam 34 from heat 3 had super high protein
LarExp=LarExp[-c(52), ]
LarExp$Log=log(LarExp$PrPerLar)
hist(LarExp$Log)

quartz()

g8=ggplot(LarExp,aes(factor(Trmt),Log,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y="Log (µg protein/larvae)")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")

all.marg2=summarySE(LarExp,measurevar="Log",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg2

pd <- position_dodge(.2)
quartz()

p4=ggplot(all.marg2,aes(x=Fam,y=Log,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(x="Family",y="Log (µg protein/larvae)",title="Larvae")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")

###Stats###
model2<- lme(Log~Trmt*Origin,data=LarExp,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model2)


# Fixed effects:  Log ~ Trmt * Origin 
# Value  Std.Error DF   t-value p-value
# (Intercept)              1.0297485 0.11592279 47  8.883055  0.0000
# TrmtHeat                -0.0483108 0.06842557 47 -0.706035  0.4837
# OriginOffshore          -0.0933719 0.16393959  8 -0.569551  0.5846
# TrmtHeat:OriginOffshore  0.0020367 0.09775556 47  0.020835  0.9835


#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model2),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model2,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model2) #0.02894155 0.6236329

##############Recruit Protein##############

RecPro=read.csv("~/Dropbox/PoritesSpawnMay2019/Protein/RecPro.csv")
str(RecPro)
RecPro$Fam=as.factor(RecPro$Fam)
RecPro$Rep=as.factor(RecPro$Rep)
hist(RecPro$PrPerSA) # Needs log transformaiton
RecPro[RecPro$PrPerSA>800,] #Fam 11 from heat 1 had super high protein
RecPro=RecPro[-c(6), ]

RecPro$Log=log(RecPro$PrPerSA)
hist(RecPro$Log)

quartz()

g9=ggplot(RecPro,aes(factor(Trmt),Log,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(Log~(µg~protein~per~cm^2)))+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(axis.title.x = element_blank())+ #remove x axis title
  theme(legend.position = "none")

all.marg3=summarySE(RecPro,measurevar="Log",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg3

pd <- position_dodge(.2)
quartz()

p6=ggplot(all.marg3,aes(x=Fam,y=Log,colour=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(Log~(µg~protein~per~cm^2)),title="Recruit")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(axis.title.x = element_blank())

###Stats###
model3<- lme(Log~Trmt*Origin,data=RecPro,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model3)

# Fixed effects:  Log ~ Trmt * Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)              4.224818 0.2204803 47 19.161878  0.0000
# TrmtHeat                -0.226389 0.1168108 47 -1.938086  0.0586
# OriginOffshore           0.200887 0.3118063  8  0.644270  0.5374
# TrmtHeat:OriginOffshore -0.169879 0.1635543 47 -1.038671  0.3043

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model3),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model3,which=1)

r.squaredGLMM(model3) # Fam: 0.08832434 0.7096161 # Fam + tank: 0.08832434 0.9464237

ggarrange(p1,p3,p5,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))
ggarrange(p2,p4,p6,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))

save(AdultPro, LarExp, RecPro, file= "AllPro.RData")
