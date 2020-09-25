rm(list=ls(all=TRUE))  # clean up left overs from any previous run
setwd("/Users/radivot/case/active/papers/TK2/R") 
dd=read.table("../data/TK2.txt",header=T)
###d=transform(dd,subset=(dT>0)&(dC==0)&(dCTP==0)&(dTTP==0)&(fig=="4")&(year==1991),select=c(dT,VdT,year,fig))
d=subset(dd,subset=(dT>0)&(dC==0)&(dCTP==0)&(dTTP==0)&(fig=="4")&(year==1991),select=c(dT,VdT,year,fig))
#d=subset(dd,subset=(dT>0)&(dC==0)&(dCTP==0)&(dTTP==0),select=c(dT,VdT,year,fig))
d
Vmic=500
S50ic=1
hillic=1
with(d, plot(dT,VdT,log="x"))
#
#MMmod<-nls(VdT~Vmax*(dT)/((dT50)+(dT)), start=list(Vmax=1000,dT50=50),trace=T,data=subset(d,year==2003))
#MMmod
hill03<-nls(VdT~Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill), 
    start=list(Vmax=Vmic,dT50=S50ic,hill=hillic),algorithm="port",data=d)
#hill03<-nls(VdT~Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill), start=list(Vmax=1000,dT50=5,hill=.5),algorithm="port",trace=T,data=subset(d,year==2003,fig=="t1"))
summary(hill03)
#confint(hill03)


# compared to JBC 2003 where Vm= 1288 +- 72, Km=13+-3, and h= .5 (no ci), 
# here we have
#      Estimate Std. Error t value Pr(>|t|)    
#Vmax 1125.6375    69.9342  16.096 8.71e-05 ***
#dT50    9.3354     1.9285   4.841  0.00840 ** 
#hill    0.8145     0.1174   6.940  0.00226 ** 

#            2.5%       97.5%
# Vmax 981.065580 1488.092496
# dT50   5.982998   26.505617
# hill   0.506972    1.202718


hill99<-nls(VdT~Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill), start=list(Vmax=1000,dT50=5,hill=.5),algorithm="port",trace=T,data=subset(d,year==1999))
hill99
confint(hill99)

# Nonlinear regression model
#   model:  VdT ~ Vmax * (dT/dT50)^hill/(1 + (dT/dT50)^hill) 
#    data:  subset(d, year == 1999) 
#     Vmax     dT50     hill 
# 863.2912  15.5436   0.5251 
#  residual sum-of-squares: 2335

#             2.5%        97.5%
# Vmax 767.9648867 1030.8279240
# dT50   9.3706969   36.2889037
# hill   0.4510934    0.6028705

# this fails to converge using the same ICs as above as well as with these
#hill91<-nls(VdT~Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill), start=list(Vmax=800,dT50=15,hill=.5),algorithm="port",trace=T,data=subset(d,year==1991))


par(mfrow=c(3,1))

lcofs<-as.list(coef(hill03));attach(lcofs)
d1=subset(d,year==2003); attach(d1)
plot(dT,VdT,type="p")
hill1=Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill)
hill2=1288*(dT/13)^.5/(1+(dT/13)^.5)
lines(dT,hill1,col="blue")
lines(dT,hill2,col="red")
VMM=coef(MMmod)["Vmax"]
KMM=coef(MMmod)["dT50"]
MM1=VMM*(dT/KMM)/(1+(dT/KMM))
lines(dT,MM1,col="green")
detach(lcofs)
detach(d1)


lcofs<-as.list(coef(hill99));attach(lcofs)
d1=subset(d,year==1999); attach(d1)
plot(dT,VdT,type="p")
hill1=Vmax*(dT/dT50)^hill/(1+(dT/dT50)^hill)
lines(dT,hill1,col="blue")
detach(lcofs)
detach(d1)

d1=subset(d,year==1991); attach(d1)
plot(dT,VdT,type="p")
hill1=700*(dT/13)^.4/(1+(dT/13)^.4) # trial and error by hand since it didn't converge
lines(dT,hill1,col="blue")
detach(d1)

