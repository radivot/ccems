# time course experiments
dT=c(0.02,0.1,0.4,2)
vol=50 # uL
dT*vol # in pM

0.001 *50 # in pm = 0.050 pmoles 

setwd("/Users/radivot/case/active/papers/TK2/data") 
stds=read.table("stds.txt",header=T)
plot(stds,log="xy")
lm1<-lm(cpm~pmol,data=stds)
lines(lm1)

d=read.table("TK2dTtimeCrs.txt",header=T)
concs=c(0.1,0.2,0.4,0.6,0.8,1,2,4,8,20,40,100,100,120,160)
Times=rep(times<-c(0,4,8,12),length(concs))
Concs=rep(concs,each=4)
d=data.frame(Times,Concs,d)
library(rgl)
plot3d(d$Times,log10(d$Concs),d$cpm)
material3d(alpha=1)
surface3d(times,log10(concs),matrix(d$cpm,ncol=length(concs)),col="violet")

d


