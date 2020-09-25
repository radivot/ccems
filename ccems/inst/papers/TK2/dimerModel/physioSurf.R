rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/ARD12.RData")
#load("C:/Users/radivot/case/active/papers/TK2/CD12.RData")
load("C:/Users/radivot/case/active/papers/TK2/CD17.RData")
#gCD
#gARD

Xf <- c(.001,.01,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50)  # fine for surfaces
Xg <- c(0.5,1,2,5,10,20,50)   # grid of experimental design points

Xf <- c(.001,.002,.005,.0075,.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100)  # fine for surfaces
Xg <- c(.001,.01,.1,1,10,100)  # fine for surfaces

fgCD=gCD
fgARD=gARD

Top=1  # T operating point
Cop=0  # C operating point
pf=data.frame(T=Top,C=Cop,t=rep(Xf,length(Xf)),c=rep(Xf,each=length(Xf)),It=1,Ic=0)
pf

fgCD$d=pf   # predicted fine mesh
fgARD$d=pf
fgCD=getEk(fgCD)
fgARD=getEk(fgARD)
fgCD
fgARD


pd=data.frame(T=Top,C=Cop,t=rep(Xg,length(Xg)),c=rep(Xg,each=length(Xg)),It=1,Ic=0)
gCD$d=pd     # predicted data (to be created)
gARD$d=pd
gCD=getEk(gCD)
gARD=getEk(gARD)

library(rgl)
open3d()

pd=gCD$d
with(pd,plot3d(log10(t),log10(c),Ek,col="violet",zlim=c(0,0.08),
        type="n",radius=.05,ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec"))
#pd=gARD$d
#with(pd,plot3d(log10(t),log10(c),Ek,col="blue",add=T,
#        type="s",radius=.05))
with(fgCD$d,surface3d(log10(Xf),log10(Xf),alpha=0.5,t(matrix(Ek[It==1],byrow=TRUE,ncol=length(Xf))),col="violet"))
with(fgARD$d,surface3d(log10(Xf),log10(Xf),alpha=0.5,t(matrix(Ek[It==1],byrow=TRUE,ncol=length(Xf))),col="blue"))

#######################  Show worse discrimination for dC ##########################


rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
source("C:/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/common.R")
load("C:/Users/radivot/case/active/papers/TK2/ARD12.RData")
load("C:/Users/radivot/case/active/papers/TK2/CD17.RData")
#load("C:/Users/radivot/case/active/papers/TK2/CD12.RData")
#gCD
#gARD

Xf <- c(.001,.01,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50)  # fine for surfaces
Xg <- c(0.5,1,2,5,10,20,50)   # grid of experimental design points

Xf <- c(.001,.002,.005,.0075,.01,.015,.02,.035,.1,.15,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,75,100)  # fine for surfaces
Xg <- c(.001,.01,.1,1,10,100)  # fine for surfaces

fgCD=gCD
fgARD=gARD



Top=0  # T operating point
Cop=1  # C operating point
pf= data.frame(T=Top,C=Cop,t=rep(Xf,length(Xf)),c=rep(Xf,each=length(Xf)),It=0,Ic=1) 
pf

fgCD$d=pf   # predicted fine mesh
fgARD$d=pf
fgCD=getEk(fgCD)
fgARD=getEk(fgARD)
fgCD
fgARD

pd=data.frame(T=Top,C=Cop,t=rep(Xg,length(Xg)),c=rep(Xg,each=length(Xg)),It=0,Ic=1) 
gCD$d=pd     # predicted data (to be created)
gARD$d=pd
gCD=getEk(gCD)
gARD=getEk(gARD)

library(rgl)
open3d()

pd=gCD$d
with(pd,plot3d(log10(t),log10(c),Ek,col="violet",
        type="n",radius=.1,ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec",zlim=c(0,.05)))
#pd=gARD$d
#with(pd,plot3d(log10(t),log10(c),Ek,col="blue",add=T,
#        type="s",radius=.1,ylab="[dCTP]",xlab="[dTTP]",zlab="1/sec",zlim=c(0,.05)))
with(fgCD$d,surface3d(log10(Xf),log10(Xf),alpha=0.5,t(matrix(Ek[Ic==1],byrow=TRUE,ncol=length(Xf))),col="violet"))
with(fgARD$d,surface3d(log10(Xf),log10(Xf),alpha=0.5,t(matrix(Ek[Ic==1],byrow=TRUE,ncol=length(Xf))),col="blue"))

