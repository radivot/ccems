rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/ARD12.RData")
gARD
Xft <- c(.0001,.0002,.0005,.00075,.001,.002,.005,.0075,.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100)  # fine for surfaces
Xfc <- c(.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100,150,200,350,500,750,1000)  # fine for surfaces
#Phyt <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,30)  # fine for surfaces
#Phyc <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20)  # fine for surfaces
Phyt <- c(3,5,7.5,10,15,20,30)  # fine for surfaces
Phyc <- c(2,3.5,5,7.5,10,15)  # fine for surfaces
fgARD=gARD
PgARD=gARD
Top=1  # T operating point
Cop=0  # C operating point
pf=data.frame(T=Top,C=Cop,t=rep(Xft,length(Xfc)),c=rep(Xfc,each=length(Xft)),It=1,Ic=0)
pf
fgARD$d=pf
fgARD=getEk(fgARD)
fgARD

Pf=data.frame(T=Top,C=Cop,t=rep(Phyt,length(Phyc)),c=rep(Phyc,each=length(Phyt)),It=1,Ic=0)
PgARD$d=Pf
PgARD=getEk(PgARD)
PgARD

library(rgl)
open3d()
pd=gARD$d
pd=subset(pd,(T==1)&(C==0))
with(pd,plot3d(log10(t+1e-4),log10(c+1e-2),Ek,xlim=c(-4,2),ylim=c(-2,3),zlim=c(0,0.082),
#        type="p",radius=.03,ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec"))
#        type="s",radius=.05,ylab="",xlab="",zlab="")) # causes Z axis problems to have data as speheres
        type="p",ylab="",xlab="",zlab=""))
with(fgARD$d,surface3d(log10(Xft),log10(Xfc),alpha=0.5,t(matrix(Ek[It==1],byrow=TRUE,ncol=length(Xfc))),col="blue"))
with(PgARD$d,surface3d(log10(Phyt),log10(Phyc),alpha=1,t(matrix(0.99*Ek[It==1],byrow=TRUE,ncol=length(Phyc))),col="red"))

#######################  Now surface for dC activity ##########################
rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/ARD12.RData")
gARD
Xft <- c(.0001,.0002,.0005,.00075,.001,.002,.005,.0075,.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100)  # fine for surfaces
Xfc <- c(.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100,150,200,350,500,750,1000)  # fine for surfaces
#Phyt <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,30)  # fine for surfaces
#Phyc <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20)  # fine for surfaces
Phyt <- c(3,5,7.5,10,15,20,30)  # fine for surfaces
Phyc <- c(2,3.5,5,7.5,10,15)  # fine for surfaces
fgARD=gARD
PgARD=gARD

Top=0  # T operating point
Cop=1  # C operating point
pf=data.frame(T=Top,C=Cop,t=rep(Xft,length(Xfc)),c=rep(Xfc,each=length(Xft)),It=0,Ic=1)
pf
fgARD$d=pf
fgARD=getEk(fgARD)
fgARD

Pf=data.frame(T=Top,C=Cop,t=rep(Phyt,length(Phyc)),c=rep(Phyc,each=length(Phyt)),It=0,Ic=1)
PgARD$d=Pf
PgARD=getEk(PgARD)
PgARD

library(rgl)
open3d()
pd=gARD$d
pd=subset(pd,(T==0)&(C==1))
with(pd,plot3d(log10(t+1e-4),log10(c+1e-2),Ek,xlim=c(-4,2),ylim=c(-2,3),zlim=c(0,0.05),
#        type="p",radius=.03,ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec"))
        type="p",ylab="",xlab="",zlab=""))
with(fgARD$d,surface3d(log10(Xft),log10(Xfc),alpha=0.5,t(matrix(Ek[Ic==1],byrow=TRUE,ncol=length(Xfc))),col="blue"))
with(PgARD$d,surface3d(log10(Phyt),log10(Phyc),alpha=1,t(matrix(0.99*Ek[Ic==1],byrow=TRUE,ncol=length(Phyc))),col="red"))

