source('summarySE.R')

#install.packages("MCMCglmm")
library(MCMCglmm)
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

##############AdultSym##############

AdultSym=read.csv("~/Dropbox/PoritesSpawnMay2019/Sym/AdtSym.csv")
str(AdultSym)
AdultSym$Fam=as.factor(AdultSym$Fam)
AdultSym$Rep=as.factor(AdultSym$Rep)
hist(AdultSym$Conc)
AdultSym[AdultSym$Conc>1000000,] #Fam 46 from control 1 had super high sym density
#AdultSym$Log=log(AdultSym$Conc)
#hist(AdultSym$Log) #still doesn't look normally distributed
AdultSym=AdultSym[-c(56), ]
#AdultSym=AdultSym[-c(23,33), ]

quartz()

g1=ggplot(AdultSym,aes(factor(Trmt),Conc,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(Cells~per~cm^2),title="Adult")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")+ #remove legend
  theme(axis.title.x = element_blank()) #remove x axis title

all.marg1=summarySE(AdultSym,measurevar="Conc",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg1

pd <- position_dodge(.2)
quartz()

p2=ggplot(all.marg1,aes(x=Fam,y=Conc,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Conc+se,ymin=Conc-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(Cells~per~cm^2),title="Adult")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) #Turn off scientific notation

###Stats###
model1<- lme(Conc~Trmt*Origin,data=AdultSym,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model1) #Sig treatment effect, heat has less
# Fixed effects:  Conc ~ Trmt * Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)              489091.8  77304.10 47  6.326855  0.0000
# TrmtHeat                -264909.5  50132.65 47 -5.284171  0.0000
# OriginOffshore            46201.8 109793.23  8  0.420808  0.6850
# TrmtHeat:OriginOffshore  -97325.7  71599.88 47 -1.359300  0.1805

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model1,which=1) 

r.squaredGLMM(model1) # Fam only: 0.378008 0.7268451 # Fam+tank: 0.3760342 0.9835694

##############LarSym##############

LarSym=read.csv("~/Dropbox/PoritesSpawnMay2019/Sym/LarSym.csv")
str(LarSym)
LarSym$Fam=as.factor(LarSym$Fam)
LarSym$Rep=as.factor(LarSym$Rep)

### Pre-exp larvae
LarDRSym=subset(LarSym,Trmt==c("DR"),select=c(Fam,Trmt,Origin,Rep,CountPerLar))
hist(LarDRSym$CountPerLar)
LarDRSym[LarDRSym$CountPerLar>16000,] #Fam19Rep2,Fam5Rep3,Fam19Rep3

#LarDRSym$Day=c(rep("Day1",6),rep("Day2",4),rep("Day1",6),rep("Day2",4),rep("Day1",6),rep("Day2",4))
#str(LarDRSym) #LarDRSym$Day is character

quartz()

ggplot(LarDRSym,aes(factor(Origin),CountPerLar))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Origin",y="Cells per larvae")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))

all.marg=summarySE(LarDRSym,measurevar="CountPerLar",groupvars=c("Fam","Origin"), na.rm=T)
all.marg

pd <- position_dodge(.2)
quartz()

ggplot(all.marg,aes(x=Fam,y=CountPerLar,shape=Origin))+
  geom_errorbar(aes(ymax=CountPerLar+se,ymin=CountPerLar-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(x="Family",y="Cells per larvae")+
  theme_bw()+ 
  scale_shape_manual(values = c(2,8))

###Stats###
model<- lme(CountPerLar~Origin,data=LarDRSym,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model) # Sig Trmt and Origin effects, almost sig interaction term where offshore bleached less under heat!!!

# Fixed effects:  CountPerLar ~ Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)    9498.004  1230.371 20  7.719625  0.0000
# OriginOffshore -147.559  1740.008  8 -0.084804  0.9345

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model) #0.0003901904 0.2872954

### Experiment
LarExpSym=subset(LarSym,Trmt==c("Control","Heat"),select=c(Fam,Trmt,Origin,Rep,CountPerLar))
hist(LarExpSym$CountPerLar)
LarExpSym[LarExpSym$CountPerLar>20000,] #Fam 20 from control 3 and Fam 46 from control 3 had super high sym density
LarExpSym=LarExpSym[-c(43,59), ] #Removed the two outliers

quartz()

g2=ggplot(LarExpSym,aes(factor(Trmt),CountPerLar,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y="Cells per larvae",title="Larvae")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank()) #remove x axis title

all.marg=summarySE(LarExpSym,measurevar="CountPerLar",groupvars=c("Origin"), na.rm=T)
all.marg 
# inshore % higher=(13146.057-8774.828)/8774.828= 0.4981555 

all.marg2=summarySE(LarExpSym,measurevar="CountPerLar",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg2

pd <- position_dodge(.2)
quartz()

p4=ggplot(all.marg2,aes(x=Fam,y=CountPerLar,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=CountPerLar+se,ymin=CountPerLar-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(x="Family",y="Cells per larvae",title="Larvae")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")

###Stats###
model2<- lme(CountPerLar~Trmt*Origin,data=LarExpSym,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model2) # Sig Trmt and Origin effects, almost sig interaction term where offshore bleached less under heat!!!

# Fixed effects:  CountPerLar ~ Trmt * Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)             14025.078  1124.553 46 12.471689  0.0000
# TrmtHeat                -1725.203   761.237 46 -2.266316  0.0282
# OriginOffshore          -5332.339  1590.358  8 -3.352917  0.0100
# TrmtHeat:OriginOffshore  2109.109  1076.552 46  1.959134  0.0562

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model2),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model2,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model2) #0.3560785 0.7010188

##############RecSym##############

RecSym=read.csv("~/Dropbox/PoritesSpawnMay2019/Sym/RecSym.csv")
str(RecSym)
RecSym$Fam=as.factor(RecSym$Fam)
RecSym$Rep=as.factor(RecSym$Rep)
hist(RecSym$CountPerSA)
RecSym$Log=log(RecSym$CountPerSA)
hist(RecSym$Log) #looks normally distributed now
RecSym[RecSym$Log<11.5,] #H39-2 had low count

quartz()

g3=ggplot(RecSym,aes(factor(Trmt),Log,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(Log~(Cells~per~cm^2)),title="Recruit")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(axis.title.x = element_blank())+ #remove x axis title
  theme(legend.position = "none")

all.marg=summarySE(RecSym,measurevar="Log",groupvars=c("Trmt","Origin"), na.rm=T)
all.marg 
# inshore % change=(12.76160-13.63783)/13.63783=-0.06424996  offshore %=(12.22461-13.79140)/13.79140=-0.1136063 Dif=-0.06424996-(-0.1136063)=0.04935634

all.marg3=summarySE(RecSym,measurevar="Log",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg3

quartz()

p6=ggplot(all.marg3,aes(x=Fam,y=Log,colour=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(Log~(Cells~per~cm^2)),title="Recruit")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(axis.title.x = element_blank())

###Stats###
model3<- lme(Log~Trmt*Origin,data=RecSym,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model3) # Sig Trmt AND interaction term. INSHORE bleached less under heat!!!

# Fixed effects:  Log ~ Trmt * Origin 
# Value Std.Error DF  t-value p-value
# (Intercept)             13.637835 0.1900703 48 71.75153  0.0000
# TrmtHeat                -0.876232 0.1708218 48 -5.12951  0.0000
# OriginOffshore           0.153560 0.2688000  8  0.57128  0.5835
# TrmtHeat:OriginOffshore -0.690551 0.2415786 48 -2.85849  0.0063

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model3),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model3,which=1)

r.squaredGLMM(model3) #Fam only: 0.5620108 0.7064494 # Fam + tank: 0.5620112 0.9579824

ggarrange(p1,p3,p5,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))
ggarrange(p2,p4,p6,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))

save(AdultSym, LarExpSym, RecSym, file= "AllSymCount.RData")
