rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/CD12.RData")
gCD
#Phyt <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,30)  # fine for surfaces
#Phyc <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20)  # fine for surfaces
Phyt <- c(3,5,7.5,10,15,20,30)  # fine for surfaces
Phyc <- c(2,3.5,5,7.5,10,15)  # fine for surfaces

PgCD=gCD
Top=0.012  # T operating point
Cop=0.005  # C operating point
Pf=data.frame(T=Top,C=Cop,t=rep(Phyt,length(Phyc)),c=rep(Phyc,each=length(Phyt)),It=1,Ic=0)
PgCD$d=Pf
PgCD=getEk(PgCD)
PgCD

library(rgl)
open3d()
#plot3d(log10(range(Phyt)),log10(range(Phyc)),c(0,0.012),type="n",ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec")
plot3d(log10(range(Phyt)),log10(range(Phyc)),c(0,0.00007),type="n",ylab="",xlab="",zlab="")
with(PgCD$d,surface3d(log10(Phyt),log10(Phyc),alpha=1,t(matrix(0.99*Ek[It==1],byrow=TRUE,ncol=length(Phyc))),col="green"))


#######################  Now surface for dC activity ##########################
rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/CD12.RData")
gCD
#Phyt <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,30)  # fine for surfaces
#Phyc <- c(.05,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20)  # fine for surfaces
Phyt <- c(3,5,7.5,10,15,20,30)  # fine for surfaces
Phyc <- c(2,3.5,5,7.5,10,15)  # fine for surfaces
PgCD=gCD
Top=0.012  # T operating point
Cop=0.005  # C operating point
Pf=data.frame(T=Top,C=Cop,t=rep(Phyt,length(Phyc)),c=rep(Phyc,each=length(Phyt)),It=0,Ic=1)
PgCD$d=Pf
PgCD=getEk(PgCD)
PgCD

library(rgl)
open3d()
#plot3d(log10(range(Phyt)),log10(range(Phyc)),c(0,0.021),type="n",ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec")
plot3d(log10(range(Phyt)),log10(range(Phyc)),c(0,0.00006),type="n",ylab="",xlab="",zlab="")
with(PgCD$d,surface3d(log10(Phyt),log10(Phyc),alpha=1,t(matrix(0.99*Ek[Ic==1],byrow=TRUE,ncol=length(Phyc))),col="green"))

