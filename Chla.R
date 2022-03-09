source('summarySE.R')

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


#####adult chl
adtchl=read.csv("~/Dropbox/PoritesSpawnMay2019/Chl/AdtChl.csv")
str(adtchl)
adtchl$Fam=as.factor(adtchl$Fam)
adtchl$Rep=as.factor(adtchl$Rep)
summary(adtchl)

adtchl<-adtchl[complete.cases(adtchl$PerArea),]
hist(adtchl$PerArea)
adtchl[adtchl$PerArea>8,] # C1-46 has high chl content
adtchl$Log=log(adtchl$PerArea) #Log transform
hist(adtchl$Log)
#adtchl=adtchl[-c(11), ]

quartz()

g4=ggplot(adtchl,aes(factor(Trmt),Log,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(Log~(µg~chla~per~cm^2)))+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")+ #remove legend
  theme(axis.title.x = element_blank()) #remove x axis title

all.marg=summarySE(adtchl,measurevar="Log",groupvars=c("Trmt"), na.rm=T)
all.marg # %change btw heat and control=(0.1581946-1.0599770)/1.0599770=-0.8507566

all.marg1=summarySE(adtchl,measurevar="Log",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg1

pd=position_dodge(.3)
quartz()

p2=ggplot(all.marg1,aes(x=Fam,y=Log,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(Log~(µg~chla~per~cm^2)),title="Adult")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())

#Mixed linear model for adult
model1 <- lme(Log~Trmt*Origin,data=adtchl,random=list(Fam=~1,Tank=~1),method="REML",na.action=na.omit)
summary(model1) #Trmt significant (heat has lower)

# Fixed effects:  Log ~ Trmt * Origin 
# Value Std.Error DF   t-value p-value
# (Intercept)              1.2242463 0.3143291 48  3.894791  0.0003
# TrmtHeat                -0.8128253 0.2103656 48 -3.863870  0.0003
# OriginOffshore          -0.3285386 0.4445285  8 -0.739072  0.4810
# TrmtHeat:OriginOffshore -0.1460511 0.2965853 48 -0.492442  0.6247

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values - Exhibit a random scatter; no apparent trendline
plot(model1,which=1) 

r.squaredGLMM(model1) #Fam: 0.256211 0.6560124 #Fam+tank: 0.2554549 0.8040357

write.table(x=all.marg1, file='adtchl.csv', append=FALSE,  quote=FALSE, row.names=FALSE, col.names=TRUE, sep=',')

#####larvae chl
larchl=read.csv("~/Dropbox/PoritesSpawnMay2019/Chl/LarChl.csv")
str(larchl)
larchl$Fam=as.factor(larchl$Fam)
larchl$Rep=as.factor(larchl$Rep)
summary(larchl)

### Pre-exp larvae
larDRchl=subset(larchl,Trmt==c("DR"),select=c(Fam,Trmt,Origin,Rep,ChlaPerLar))
hist(larDRchl$ChlaPerLar)
larDRchl$Log=log(larDRchl$ChlaPerLar)
hist(larDRchl$Log)

quartz()

ggplot(larDRchl,aes(factor(Origin),Log))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Origin",y="Log (µg Chla per larvae)")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))

all.marg=summarySE(larDRchl,measurevar="Log",groupvars=c("Fam","Origin"), na.rm=T)
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
  labs(x="Family",y="Log (µg Chla per larvae)")+
  theme_bw()+ 
  scale_shape_manual(values = c(2,8))

###Stats###
model<- lme(Log~Origin,data=larDRchl,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model) 

# Fixed effects:  Log ~ Origin 
# Value Std.Error DF    t-value p-value
# (Intercept)    -4.267711 0.1671337 16 -25.534709  0.0000
# OriginOffshore -0.224756 0.2339788  8  -0.960584  0.3649

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model) #0.0731857 0.7574807 #Random effect is HUGE!!!


###Experiment
larexp=subset(larchl,Trmt==c("Control","Heat"),select=c(Fam,Trmt,Origin,Rep,ChlaPerLar))
hist(larexp$ChlaPerLar)
larexp[larexp$ChlaPerLar>0.025,] #C3-20 has too high chl
larexp=larexp[-c(43), ]

g5=ggplot(larexp,aes(factor(Trmt),ChlaPerLar,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(x="Treatment",y="µg Chla per larvae")+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank()) #remove x axis title

all.marg=summarySE(larexp,measurevar="ChlaPerLar",groupvars=c("Origin"), na.rm=T)
all.marg # %change btw inshore and offshore=(0.01838005-0.01230576)/0.01230576=0.4936136

all.marg2=summarySE(larexp,measurevar="ChlaPerLar",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg2

pd=position_dodge(0.3)
quartz()

p4=ggplot(all.marg2,aes(x=Fam,y=ChlaPerLar,color=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=ChlaPerLar+se,ymin=ChlaPerLar-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(x="Family",y="µg Chla per larvae",title="Larvae")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(legend.position = "none")

#Mixed linear model for lar

model2 <- lme(ChlaPerLar~Trmt*Origin,data=larexp,random= ~1|Fam,method="REML",na.action=na.omit)
summary(model2) #Origin significant (offshore has lower conc)

# Fixed effects:  ChlaPerLar ~ Trmt * Origin 
# Value   Std.Error DF   t-value p-value
# (Intercept)              0.019331980 0.001791721 47 10.789612  0.0000
# TrmtHeat                -0.001879043 0.001203583 47 -1.561208  0.1252
# OriginOffshore          -0.006566281 0.002522641  8 -2.602940  0.0315
# TrmtHeat:OriginOffshore  0.000959175 0.001685350 47  0.569125  0.5720

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model2),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model2,which=1)

#How good is this model overall? the proportion of variance explained by the fixed factors alone, and by both fixed and random factors.
r.squaredGLMM(model2) #0.3027784 0.6797567

write.table(x=all.marg1, file='SummarySE_larchl.csv', append=FALSE,  quote=FALSE, row.names=FALSE, col.names=TRUE, sep=',')

#####recruit chl
recchl=read.csv("~/Dropbox/PoritesSpawnMay2019/Chl/RecChl.csv")
recchl=recchl[complete.cases(recchl$ChlaPerSA),]
str(recchl)
recchl$Fam=as.factor(recchl$Fam)
recchl$Rep=as.factor(recchl$Rep)
hist(recchl$ChlaPerSA) #Needs log transformation
recchl[recchl$ChlaPerSA>8,] # C1-38 and C1-39 have high chl content

recchl$Log=log(recchl$ChlaPerSA)
#recchl[recchl$Log>2,]
#recchl=recchl[-c(13,30), ]
hist(recchl$Log)

quartz()

g6=ggplot(recchl,aes(factor(Trmt),Log,fill=Origin))+
  geom_boxplot(outlier.shape=16, notch=FALSE)+
  #geom_point(aes(colour=factor(Fam)),position=pd)+
  #geom_line(aes(group=factor(Fam),colour=factor(Fam)),position=pd)+
  labs(colour="Origin")+
  #scale_colour_manual(values=c(point="gray60"))+
  labs(y=expression(Log~(µg~chla~per~cm^2)))+
  #theme(axis.text.x=element_text(size=20))+
  theme_classic()+
  scale_fill_manual(values=c("red3","royalblue1"))+
  theme(axis.text=element_text(size=11))+
  theme(plot.subtitle=element_text(hjust=0.5))+
  theme(axis.title.x = element_blank())

all.marg=summarySE(recchl,measurevar="ChlaPerSA",groupvars=c("Trmt"), na.rm=T)
all.marg # %change btw heat and control=(-0.3803708-0.3428222)/0.3428222=-2.109528

all.marg3=summarySE(recchl,measurevar="Log",groupvars=c("Fam","Trmt","Origin"), na.rm=T)
all.marg3

quartz()

p6=ggplot(all.marg3,aes(x=Fam,y=Log,colour=Trmt,shape=Origin))+
  geom_errorbar(aes(ymax=Log+se,ymin=Log-se),lwd=0.3,width=1,position=pd)+
  #geom_line(aes(group=Family,linecol=Family),position=pd)+
  geom_point(aes(color=Trmt,shape=Origin),position=pd,size=2.5)+
  #scale_shape_identity()+
  #geom_bar(stat="identity")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=20))+
  labs(y=expression(Log~(µg~chla~per~cm^2)),title="Recruit")+
  theme_bw()+ 
  scale_color_manual(values=c("turquoise","orange"))+
  scale_shape_manual(values = c(2,8))+
  theme(axis.title.x = element_blank())

model3<- lme(Log~Trmt*Origin,data=recchl,random=list(Fam=~1),method="REML",na.action=na.omit)
summary(model3) #Almost sig trmt

# Fixed effects:  Log ~ Trmt * Origin 
# Value Std.Error DF    t-value p-value
# (Intercept)              0.3785704 0.2105864 45  1.7976966  0.0789
# TrmtHeat                -0.5962339 0.3030856 45 -1.9672132  0.0553
# OriginOffshore           0.0561762 0.3030856  8  0.1853476  0.8576
# TrmtHeat:OriginOffshore -0.3815907 0.4323070 45 -0.8826846  0.3821

#Diagnostics to determine residual is normally distributed- all good
qqPlot(residuals(model3),xlab="Theoretical Quantiles",ylab="Observed Quantiles")

#Check homoscedascitity by plotting residuals against fitted values
plot(model3,which=1)

r.squaredGLMM(model3) #0.2043329 0.2043329 # Fam+Tank 0.2043329 0.999979

ggarrange(p1,p3,p5,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))
ggarrange(p2,p4,p6,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))

save(adtchl, larexp, recchl, file= "AllChla.RData")
