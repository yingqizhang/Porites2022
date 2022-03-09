setwd("~/Dropbox/PoritesSpawnMay2019/Recruit_photos/")
source('~/Dropbox/PoritesSpawnYingqi/summarySE.R')

#install .packages(MCMCglmm)
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
#library(lmmfit)
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
library(ggplot2)
#install.packages("egg")
library(egg)
#install.packages("MuMIn")
library(MuMIn)

cs=read.csv("~/Dropbox/PoritesSpawnMay2019/Recruit_photos/ColorScore.csv")
str(cs)
cs$ID=as.factor(cs$ID)
cs$Fam=as.factor(cs$Fam)

cs$change=cs$T28-cs$T0
cs$percentchange=(cs$T28-cs$T0)/cs$T0

csCC=cs[complete.cases(cs$change),]
hist(csCC$percentchange)

recruitcs=subset(csCC,Type=="R",select=c(Fam,Trmt,Origin,Tank,change,percentchange))
hist(recruitcs$change)
adultcs=subset(csCC,Type=="A",select=c(Fam,Trmt,Origin,Tank,change,percentchange))
hist(adultcs$change)

#adult cs
quartz()

p1=ggplot(adultcs,aes(factor(Trmt),change,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y="change in color score",title="Adult")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none") #remove legend

all.marg=summarySE(adultcs,measurevar="change",groupvars=c("Trmt"), na.rm=T)
all.marg #heat reduced=(-1.2136885-(-0.3304776)/(-0.3304776))=-2.213688

all.marg1=summarySE(adultcs,measurevar="change",groupvars=c("Fam","Trmt","Origin"), na.rm=T) #change to percentchange

pd=position_dodge(.3)
quartz()

p2=ggplot(all.marg1,aes(x=Fam,y=change,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=change+se,ymin=change-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y="change in color score",title="Adult")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")

#Mixed linear model for adult
model1 <- lme(change~Trmt*Origin,data=adultcs,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model1) #Trmt significant (heat has lower)

# Fixed effects:  change ~ Trmt * Origin 
# Value Std.Error DF    t-value p-value
# (Intercept)             -0.4436952 0.3663555 48 -1.2111057  0.2318
# Trmtheat                -0.7427771 0.3173272 48 -2.3407296  0.0234
# Originoffshore           0.2264353 0.5181049  8  0.4370452  0.6736
# Trmtheat:Originoffshore -0.3241890 0.4486980 48 -0.7225105  0.4735

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model1,which=1) 

r.squaredGLMM(model1) #Fam: 0.1508832 0.45829 # Fam+Tank: 0.1555887 0.9798807

#recruit cs
quartz()

p3=ggplot(recruitcs,aes(factor(Trmt),change,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y="change in color score",title="Recruit")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))

all.marg=summarySE(recruitcs,measurevar="change",groupvars=c("Trmt"), na.rm=T)
all.marg #heat reduced=(-1.2091939-(-0.7337418)/(-0.7337418))=-2.209194

all.marg2=summarySE(recruitcs,measurevar="change",groupvars=c("Fam","Trmt","Origin"), na.rm=T) #change to percentchange

p4=ggplot(all.marg2,aes(x=Fam,y=change,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=change+se,ymin=change-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y="change in color score",title="Recruit")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))

#Mixed linear model for recruit
model2 <- lme(change~Trmt*Origin,data=recruitcs,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model2) #Trmt significant (heat has lower)

# Fixed effects:  change ~ Trmt * Origin 
# Value Std.Error  DF   t-value p-value
# (Intercept)             -0.7782541 0.1472610 451 -5.284864  0.0000
# Trmtheat                -0.3559805 0.1713889  48 -2.077033  0.0432
# Originoffshore          -0.0359505 0.2083874   8 -0.172518  0.8673
# Trmtheat:Originoffshore -0.1351730 0.2412434  48 -0.560318  0.5779

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model2),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model2,which=1) 

r.squaredGLMM(model2) #Fam:0.08536927 0.1685807 #Fam + tank: 0.07519254 0.3705669

ggarrange(p1,p3,ncol=2,labels=c('a','b'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))
ggarrange(p2,p4,ncol=2,labels=c('a','b'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))

###Recruit size
size=read.csv("~/Dropbox/PoritesSpawnMay2019/Recruit_photos/RecSize.csv")

str(size)
size$ID=as.factor(size$ID)
size$Fam=as.factor(size$Fam)

size$change=size$T28-size$T0.cm2.
size$percentchange=(size$T28-size$T0.cm2.)/size$T0.cm2.

sizeCC=size[complete.cases(size$change),]
hist(sizeCC$change)

recsize=subset(sizeCC,Type=="R",select=c(Fam,Trmt,Origin,Tank,change,percentchange))
hist(recsize$change)

#adult cs
quartz()

ggplot(recsize,aes(factor(Trmt),change,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y=expression("Change in size"~(cm^2)))+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))

all.marg=summarySE(recsize,measurevar="change",groupvars=c("Trmt"), na.rm=T)
all.marg #heat reduced=(-0.009213178-(-0.016043651)/(-0.016043651))=-1.009213

all.marg1=summarySE(recsize,measurevar="change",groupvars=c("Fam","Trmt","Origin"), na.rm=T) #change to percentchange

pd=position_dodge(.3)
quartz()

ggplot(all.marg1,aes(x=Fam,y=change,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=change+se,ymin=change-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression("Change in size"~(cm^2)))+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))

#Mixed linear model
model1 <- lme(change~Trmt*Origin,data=recsize,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model1) #Trmt significant (heat has lower)

# Fixed effects:  change ~ Trmt * Origin 
# Value   Std.Error  DF   t-value p-value
# (Intercept)             -0.016184829 0.001881985 450 -8.599874  0.0000
# Trmtheat                 0.007353238 0.002408640  48  3.052858  0.0037
# Originoffshore          -0.000281593 0.002667131   8 -0.105579  0.9185
# Trmtheat:Originoffshore -0.000418397 0.003391694  48 -0.123359  0.9023

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model1,which=1) 

r.squaredGLMM(model1) #Fam: 0.07858589 0.1385583  # Fam + tank: 0.08685059 0.3130034

