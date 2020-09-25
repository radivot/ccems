rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,3)]=c("T","C","X")
dt=subset(d,(It==1)&(T==100))
dc=subset(d,(Ic==1)&(C==100))
d=rbind(dt,dc)
d

Xg <- c(2,3.5,5,7.5,10,15,20,35,50,100,150,200,500,750,1000,2000,3500)
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(1,2),mar=c(4.2,4,.2,1)+.1,oma=c(0,0,.6,0))
with(d[(d$It>0),],plot(X,k,ylab="dT Kinase activity (1/sec)",xlab="[ATP] uM  ([dT]=100 uM)",log="x"),
    ylim=range(d$k),xlim=range(Xg))
with(d[(d$Ic>0),],plot(X,k,ylab="dC Kinase activity (1/sec)",xlab="[ATP] uM  ([dC]=100 uM)",log="x"),
    ylim=range(d$k),xlim=range(Xg))


