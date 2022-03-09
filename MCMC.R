source('summarySE.R')

library(lme4)
library(ggplot2) 
library(MCMCglmm)
library(nlme)
library(AICcmodavg)
library(MASS)
library(VGAM)
library(ordinal)
#install.packages("ordinal")
library(lattice)
library(car)
library(MCMC.qpcr)
#install.packages("MCMC.qpcr")
library(adegenet)
#install.packages("adegenet")
library(SuppDists)
library(plyr)
#install.packages("SuppDists")

se=function(x) sd(x)/sqrt(length(x))

######Sym count

load(file = "~/Dropbox/PoritesSpawnMay2019/Sym/AllSymCount.RData") #Recruit count was log transformed

##Adult
prior1=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002),G2 = list(V = 1,nu=0.002)))

#prior1=list(R =  list(V = 1, nu = 0.002),G = list(G1 = list(V = 1,nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior2=list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu=0.002)))

m1=MCMCglmm(Conc~Origin*Trmt,random=~Fam+Tank,data=AdultSym,prior=prior1,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m1)

autocorr.diag(m1$Sol) #for fixed effects
autocorr.diag(m1$VCV) #for random effects
effectiveSize(m1$Sol)
effectiveSize(m1$VCV)

summary(m1)

# Location effects: Conc ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC  
# (Intercept)                267103    74406   445821   3228.0 0.0196 *
#   OriginOffshore             100132   -54670   261858   4248.6 0.2129  
# TrmtHeat                  -193839  -310298   -20750    592.2 0.0644 .
# OriginOffshore:TrmtHeat   -102361  -212214     8263   4500.0 0.0711 .

vcv=colnames(m1$VCV) #"Fam" "Tank" "units" 
H1=(m1$VCV[,1])/(m1$VCV[,1]+m1$VCV[,2]+m1$VCV[,3]) 
mean(H1) #0.6530599
HPDinterval(H1) #0.2675134 0.9624618
mean(m1$VCV[,1]) #63242108496 #Vg
mean(m1$VCV[,2]+m1$VCV[,3]) #26897226217 #Ve

##DRlarvae
#save(LarDRSym, larDRchl, LarDRPro, file= "DRlarvae.RData")
load(file = "~/Dropbox/PoritesSpawnMay2019/MCMC/DRlarvae.RData")

m0=MCMCglmm(CountPerLar~Origin,random=~Fam,data=LarDRSym,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m0)

autocorr.diag(m0$Sol) #for fixed effects
autocorr.diag(m0$VCV) #for random effects
effectiveSize(m0$Sol)
effectiveSize(m0$VCV)

summary(m0)
# Location effects: CountPerLar ~ Origin 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)       9487.4   7294.9  11596.9     4500 <2e-04 ***
#   OriginOffshore    -119.3  -3198.2   2823.6     4500  0.936 

vcv=colnames(m0$VCV) #"Fam"   "units"
H0=(m0$VCV[,1])/(m0$VCV[,1]+m0$VCV[,2]) #Fam/sum(Fam+units)
mean(H0) #0.06891424
HPDinterval(H0) #1.298519e-11 0.4252373
mean(m0$VCV[,1]) #1082219 #Vg
mean(m0$VCV[,2]) #14180960 #Ve

##Larvae
m2=MCMCglmm(CountPerLar~Origin*Trmt,random=~Fam,data=LarExpSym,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m2)

autocorr.diag(m2$Sol) #for fixed effects
autocorr.diag(m2$VCV) #for random effects
effectiveSize(m2$Sol)
effectiveSize(m2$VCV)

summary(m2)
# Location effects: CountPerLar ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)              14023.19 11553.59 16500.41     4500 <2e-04 ***
#   OriginOffshore           -5361.88 -8735.00 -1914.11     3946 0.0080 ** 
#   TrmtHeat                 -1727.05 -3195.44  -147.10     4500 0.0298 *  
#   OriginOffshore:TrmtHeat   2101.18    65.01  4427.63     4500 0.0582 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vcv=colnames(m2$VCV) #"Fam"   "Rep"
H2=(m2$VCV[,1])/(m2$VCV[,1]+m2$VCV[,2]) #Fam/sum(Fam+units)
mean(H2) #0.5369144
HPDinterval(H2) #0.2549022 0.8384471
mean(m2$VCV[,1]) #6303971 #Vg
mean(m2$VCV[,2]) #4386180 #Ve

##Recruit
m3=MCMCglmm(Log~Origin*Trmt,random=~Fam+Tank,data=RecSym,prior=prior1,nitt=100000,thin=20,burnin=10000)

summary(m3)
# Location effects: Log ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp    pMCMC    
# (Intercept)               13.6336  13.1772  14.0814     4500  < 2e-04 ***
#   OriginOffshore             0.1554  -0.4302   0.7160     5200 0.578222    
# TrmtHeat                  -0.8749  -1.2866  -0.4348     5098 0.000444 ***
#   OriginOffshore:TrmtHeat   -0.6931  -1.1636  -0.1998     4500 0.006222 ** 

vcv=colnames(m3$VCV) #"Fam" "Tank"  "units"
H3=(m3$VCV[,1])/(m3$VCV[,1]+m3$VCV[,2]+m3$VCV[,3])
mean(H3) #0.315865
HPDinterval(H3) #0.0008881933 0.6102024
mean(m3$VCV[,1]) #0.1359263 #Vg
mean(m3$VCV[,2]+m3$VCV[,3]) #0.2536101 #Ve

AdtSym=c(0.6530599,0.2675134,0.9624618)
LarDRSym=c(0.06891424,1.298519e-11,0.4252373)
LarSym=c(0.5369144,0.2549022,0.8384471)
RecSym=c(0.315865,0.0008881933,0.6102024)

Sym=rbind(AdtSym,LarDRSym,LarSym,RecSym)
colnames(Sym)=c("mean","lwr","upper")
Sym=data.frame(Sym)
Sym$type=c("Adult","LarvaeT0","LarvaeT1","Recruit")
Sym$trait=c("Symbiont","Symbiont","Symbiont","Symbiont")

p1=ggplot(Sym, aes(x=factor(type,level=c("Adult","LarvaeT0","LarvaeT1","Recruit")), y=mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=lwr, ymax=upper),
                size=.3,    # Thinner lines
                width=.2) +
  theme_bw()+
  labs(y=expression(H^2),title="Symbiont")+
  theme(axis.title.x = element_blank()) #remove x axis title

######Chl conc

load(file = "~/Dropbox/PoritesSpawnMay2019/Chl/AllChla.RData") #Adult and Recruit Chl were log transformed

##Adult
m1=MCMCglmm(Log~Origin*Trmt,random=~Fam+Tank,data=adtchl,prior=prior1,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m1)

autocorr.diag(m1$Sol) #for fixed effects
autocorr.diag(m1$VCV) #for random effects
effectiveSize(m1$Sol)
effectiveSize(m1$VCV)

summary(m1)

# Location effects: Log ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp   pMCMC   
# (Intercept)                1.2262   0.5145   1.9554     4500 0.00489 **
#   OriginOffshore            -0.3319  -1.3292   0.6432     5394 0.47289   
# TrmtHeat                  -0.8143  -1.3385  -0.2978     4500 0.00756 **
#   OriginOffshore:TrmtHeat   -0.1455  -0.7039   0.4674     4500 0.61956   

vcv=colnames(m1$VCV) #"Fam" "Tank"  "units"
H1=(m1$VCV[,1])/(m1$VCV[,1]+m1$VCV[,2]+m1$VCV[,3])
mean(H1) #0.5332771
HPDinterval(H1) #0.2654002 0.8378967
mean(m1$VCV[,1]) #0.5027665 #Vg
mean(m1$VCV[,2]+m1$VCV[,3]) #0.379871 #Ve

##DRlarvae
m0=MCMCglmm(Log~Origin,random=~Fam,data=larDRchl,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m0)

autocorr.diag(m0$Sol) #for fixed effects
autocorr.diag(m0$VCV) #for random effects
effectiveSize(m0$Sol)
effectiveSize(m0$VCV)

summary(m0)
# Location effects: Log ~ Origin 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)      -4.2657  -4.6221  -3.8854     4500 <2e-04 ***
#   OriginOffshore   -0.2264  -0.7185   0.3263     4500   0.36    

vcv=colnames(m0$VCV) #"Fam"   "units"
H0=(m0$VCV[,1])/(m0$VCV[,1]+m0$VCV[,2]) #Fam/sum(Fam+units)
mean(H0) #0.7137377
HPDinterval(H0) #0.4282521 0.9564523
mean(m0$VCV[,1]) #0.1562911 #Vg
mean(m0$VCV[,2]) #0.05012536 #Ve

##Larvae
m2=MCMCglmm(ChlaPerLar~Origin*Trmt,random=~Fam,data=larexp,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m2)

autocorr.diag(m2$Sol) #for fixed effects
autocorr.diag(m2$VCV) #for random effects
effectiveSize(m2$Sol)
effectiveSize(m2$VCV)

summary(m2)
# Location effects: ChlaPerLar ~ Origin * Trmt 
# 
# post.mean   l-95% CI   u-95% CI eff.samp  pMCMC  
# (Intercept)              0.0192758  0.0017515  0.0356596     4500 0.0316 *
#   OriginOffshore          -0.0063831 -0.0321111  0.0164764     4500 0.5693  
# TrmtHeat                -0.0017807 -0.0073691  0.0035647     4500 0.5169  
# OriginOffshore:TrmtHeat  0.0008224 -0.0067025  0.0084263     3858 0.8169  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vcv=colnames(m2$VCV) #"Fam"   "units"
H2=(m2$VCV[,1])/(m2$VCV[,1]+m2$VCV[,2]) #Fam/sum(Fam+units)
mean(H2) #0.8372514
HPDinterval(H2) #0.6962491 0.9576276
mean(m2$VCV[,1]) #0.0003640134 #Vg
mean(m2$VCV[,2]) #5.488132e-05 #Ve

##Recruit
m3=MCMCglmm(Log~Origin*Trmt,random=~Fam+Tank,data=recchl,prior=prior1,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m3) #density of Fam looks wonky, might be due to fewer entries

summary(m3)
# Location effects: Log ~ Origin * Trmt 
# 
# post.mean  l-95% CI  u-95% CI eff.samp pMCMC  
# (Intercept)              0.363600 -0.874049  1.826119     4500 0.495  
# OriginOffshore           0.112985 -0.348026  0.581920     4500 0.604  
# TrmtHeat                -0.474098 -2.439067  1.407633     4089 0.528  
# OriginOffshore:TrmtHeat -0.471759 -0.927192  0.006245     4902 0.056 .

vcv=colnames(m3$VCV) #"Fam" "Tank"  "units"
H3=(m3$VCV[,1])/(m3$VCV[,1]+m3$VCV[,2]+m3$VCV[,3])
mean(H3) #0.07189849
HPDinterval(H3) #4.027399e-05 0.222543
mean(m3$VCV[,1]) #0.07711243 #Vg
mean(m3$VCV[,2]+m3$VCV[,3]) #1.782732 #Ve

AdtChl=c(0.5332771,0.2654002,0.8378967)
LarDRChl=c(0.7137377,0.4282521,0.9564523)
LarChl=c(0.8372514,0.6962491,0.9576276)
RecChl=c(0.07189849,4.027399e-05,0.222543)

Chl=rbind(AdtChl,LarDRChl,LarChl,RecChl)
colnames(Chl)=c("mean","lwr","upper")
Chl=data.frame(Chl)
Chl$type=c("Adult","LarvaeT0","LarvaeT1","Recruit")
Chl$trait=c("Chla","Chla","Chla","Chla")


p2=ggplot(Chl, aes(x=factor(type,level=c("Adult","LarvaeT0","LarvaeT1","Recruit")), y=mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=lwr, ymax=upper),
                size=.3,    # Thinner lines
                width=.2) +
  theme_bw()+
  labs(x="Life stage",title="Chl a")+
  theme(axis.title.y = element_blank()) #remove y axis title
 
######Protein

load(file = "~/Dropbox/PoritesSpawnMay2019/Protein/AllPro.RData") #Lar and Recruit protein were log transformed

##Adult

m1=MCMCglmm(PerArea~Origin*Trmt,random=~Fam+Tank,data=AdultPro,prior=prior1,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m1)

autocorr.diag(m1$Sol) #for fixed effects
autocorr.diag(m1$VCV) #for random effects
effectiveSize(m1$Sol)
effectiveSize(m1$VCV)

summary(m1)

# Location effects: PerArea ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)               118.305   94.975  141.287     4500 <2e-04 ***
#   OriginOffshore             10.848  -18.837   43.369     4623 0.4498    
# TrmtHeat                   -1.275  -23.565   20.675     4500 0.8884    
# OriginOffshore:TrmtHeat   -29.697  -57.751   -1.730     4500 0.0364 * 

vcv=colnames(m1$VCV) #"Fam"  "Tank" "units"
H1=(m1$VCV[,1])/(m1$VCV[,1]+m1$VCV[,2]+m1$VCV[,3])
mean(H1) #0.2620026
HPDinterval(H1) #2.542921e-07 0.5792513
mean(m1$VCV[,1]) #341.1442 #Vg
mean(m1$VCV[,2]+m1$VCV[,3]) #803.7168 #Ve

##DRlarvae
m0=MCMCglmm(Log~Origin,random=~Fam,data=LarDRPro,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m0)

autocorr.diag(m0$Sol) #for fixed effects
autocorr.diag(m0$VCV) #for random effects
effectiveSize(m0$Sol)
effectiveSize(m0$VCV)

summary(m0)
# Location effects: Log ~ Origin 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)      1.20150  0.96484  1.42665     4500 <2e-04 ***
#   OriginOffshore  -0.08805 -0.39009  0.25642     4500  0.537     

vcv=colnames(m0$VCV) #"Fam"   "units"
H0=(m0$VCV[,1])/(m0$VCV[,1]+m0$VCV[,2]) #Fam/sum(Fam+units)
mean(H0) #0.3701773
HPDinterval(H0) #0.004409066 0.7274645
mean(m0$VCV[,1]) #0.0434072 #Vg
mean(m0$VCV[,2]) #0.05834261 #Ve

##Larvae
m2=MCMCglmm(Log~Origin*Trmt,random=~Fam,data=LarExp,prior=prior2,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m2)

autocorr.diag(m2$Sol) #for fixed effects
autocorr.diag(m2$VCV) #for random effects
effectiveSize(m2$Sol)
effectiveSize(m2$VCV)

summary(m2)
# Location effects: Log ~ Origin * Trmt 
# 
# post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)              1.033579  0.773665  1.287990     4865 <2e-04 ***
#   OriginOffshore          -0.095163 -0.469947  0.262099     4767  0.564    
# TrmtHeat                -0.048714 -0.181851  0.093807     4753  0.489    
# OriginOffshore:TrmtHeat  0.003115 -0.195004  0.201224     4711  0.972    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vcv=colnames(m2$VCV) #"Fam"   "Rep"
H2=(m2$VCV[,1])/(m2$VCV[,1]+m2$VCV[,2]) #Fam/sum(Fam+units)
mean(H2) #0.6103649
HPDinterval(H2) #0.3327874 0.8654167
mean(m2$VCV[,1]) #0.07154907 #Vg
mean(m2$VCV[,2]) #0.03682562 #Ve

##Recruit
m3=MCMCglmm(Log~Origin*Trmt,random=~Fam+Tank,data=RecPro,prior=prior1,nitt=100000,thin=20,burnin=10000)

quartz()
plot(m3) 

summary(m3)
# Location effects: Log ~ Origin * Trmt 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)               4.22248  3.69440  4.72181     4500 <2e-04 ***
#   OriginOffshore            0.20196 -0.54986  0.87064     4500  0.518    
# TrmtHeat                 -0.22527 -0.53692  0.09667     4283  0.135    
# OriginOffshore:TrmtHeat  -0.16808 -0.46264  0.16740     4500  0.295    

vcv=colnames(m3$VCV) 
H3=(m3$VCV[,1])/(m3$VCV[,1]+m3$VCV[,2]+m3$VCV[,3]) 
mean(H3) #0.666978
HPDinterval(H3) #0.4235869 0.9153994
mean(m3$VCV[,1]) #0.2780672 #Vg
mean(m3$VCV[,2]+m3$VCV[,3]) #0.1157802 #Ve

AdtPro=c(0.2620026,2.542921e-07,0.5792513)
LarDRPro=c(0.3701773,0.004409066,0.7274645)
LarPro=c(0.6103649,0.3327874,0.8654167)
RecPro=c(0.666978,0.4235869,0.9153994)

Pro=rbind(AdtPro,LarDRPro,LarPro,RecPro)
colnames(Pro)=c("mean","lwr","upper")
Pro=data.frame(Pro)
Pro$type=c("Adult","LarvaeT0","LarvaeT1","Recruit")
Pro$trait=c("Protein","Protein","Protein","Protein")

p3=ggplot(Pro, aes(x=factor(type,level=c("Adult","LarvaeT0","LarvaeT1","Recruit")), y=mean)) + 
  geom_bar(stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=lwr, ymax=upper),
                size=.3,    # Thinner lines
                width=.2) +
  theme_bw()+
  labs(title="Protein")+
  theme(axis.title.x = element_blank())+ #remove x axis title
  theme(axis.title.y = element_blank()) #remove y axis title
  

quartz()

ggarrange(p1,p2,p3,ncol=3,labels=c('a','b','c'),label.args=list(gp=grid::gpar(font=1,cex=1.2)))

###Color score
cs=read.csv("~/Dropbox/PoritesSpawnMay2019/Recruit_photos/ColorScore.csv")
str(cs)
cs$ID=as.factor(cs$ID)
cs$Fam=as.factor(cs$Fam)

cs$change=cs$T28-cs$T0
csCC=cs[complete.cases(cs$change),]
histogram(csCC$change)

#Analyze adult and recruit separately
recruitcs=subset(csCC,Type=="R",select=c(Fam,Trmt,Origin,Tank,change))
adultcs=subset(csCC,Type=="A",select=c(Fam,Trmt,Origin,Tank,change))

mR=MCMCglmm(change~Origin*Trmt,random=~Fam+Tank,data=recruitcs,prior=prior1,nitt=100000,thin=20,burnin=10000)
summary(mR)
HR=(mR$VCV[,1])/(mR$VCV[,1]+mR$VCV[,2]+mR$VCV[,3])
mean(HR) #0.1064759
HPDinterval(HR) #0.01426336 0.2512983
mean(mR$VCV[,1]) #0.07569104
mean(mR$VCV[,2]+mR$VCV[,3]) #0.5987496

mA=MCMCglmm(change~Origin*Trmt,random=~Fam+Tank,data=adultcs,prior=prior1,nitt=100000,thin=20,burnin=10000)
summary(mA)
HA=(mA$VCV[,1])/(mA$VCV[,1]+mA$VCV[,2]+mA$VCV[,3])
mean(HA) #0.3537751
HPDinterval(HA) #0.03636969 0.666462
mean(mA$VCV[,1]) #0.5468923
mean(mA$VCV[,2]+mA$VCV[,3]) #0.8698956

#Rec size
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

m1=MCMCglmm(change~Origin*Trmt,random=~Fam+Tank,data=recsize,prior=prior1,nitt=100000,thin=20,burnin=10000)
summary(m1)

# Location effects: change ~ Origin * Trmt 
# 
# post.mean   l-95% CI   u-95% CI eff.samp pMCMC
# (Intercept)             -0.0158909 -0.0550831  0.0254131     4500 0.359
# Originoffshore          -0.0009650 -0.0248702  0.0220583     4500 0.927
# Trmtheat                 0.0066026 -0.0440547  0.0585801     4500 0.756
# Originoffshore:Trmtheat  0.0006448 -0.0030023  0.0042422     4500 0.738

vcv=colnames(m1$VCV)
H1=(m1$VCV[,1])/(m1$VCV[,1]+m1$VCV[,2]+m1$VCV[,3])
mean(H1) #0.2946306
HPDinterval(H1) #0.03270691 0.5992974
