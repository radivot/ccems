rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(year==2003)&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
dtm=subset(d,(It==1)&(C==0))
dcm=subset(d,(Ic==1)&(T==0))
d=rbind(dtm,dcm)
if (DNA=="wt") {
  tinh=list(list(t=0,c=c(0,1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0))
  cinh=list(list(t=0,c=c(0,1,5)),list(t=1,c=0),list(t=5,c=0)) } else { 
  tinh=list(list(t=0,c=c(0,1,10)),list(t=1,c=0),list(t=10,c=0))
  cinh=list(list(t=0,c=c(0,2,10)),list(t=2,c=0),list(t=10,c=0)) 
}
Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)

if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(1,2),mar=c(4.2,4,.2,1)+.1,oma=c(0,0,1.2,0))

plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
legs=NULL
for (i in 1:length(tinh)) 
{
for (j in 1:length(tinh[[i]]$c)) {
  with(d[(d$It>0)&(d$t==tinh[[i]]$t)&(d$c==tinh[[i]]$c[j]),],points(T,k,col=i,pch=j))
  legs=rbind(legs,list(s=paste("[dTTP] =",format(tinh[[i]]$t,width=3),";  [dCTP] =",tinh[[i]]$c[j],sep=""),c=i,p=j))
  }
}
legs=as.data.frame(legs)
legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")

plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase activity (1/sec)",xlab="[dC] uM",log="x")
legs=NULL
for (i in 1:length(cinh)) 
{
for (j in 1:length(cinh[[i]]$c)) {
  with(d[(d$Ic>0)&(d$t==cinh[[i]]$t)&(d$c==cinh[[i]]$c[j]),],points(C,k,col=i,pch=j))
  legs=rbind(legs,list(s=paste("[dTTP] =",format(cinh[[i]]$t,width=3),";  [dCTP] =",cinh[[i]]$c[j],sep=""),c=i,p=j))
  }
}
legs=as.data.frame(legs)
legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
title(DNA,outer=TRUE,line=0.2,cex.main=1)


