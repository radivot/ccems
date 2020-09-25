rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 12, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(2,6),mar=c(4,4,2,1)+.1)

plotCol<-function(MM2,dd,xlabS) {
plot(dd$SF,dd$k,xlab=paste("Total [",xlabS,"]",sep=""),log="xy", ylab="k (1/sec) ", 
    main=paste(d[1,"frstAut"],d[1,"year"]))
lgx=log(dd$SF)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
lines(fineX,predict(MM2,list(SF=fineX)),col="black",lwd=1)
plot(MM2$m$fitted(),MM2$m$resid(),xlab="Fitted Value",
    ylab="Relative Residual",mar=c(2,2,0,1)+.1)
}

n=1:3
(1+(1-1/n)^0.5)/(1-(1-1/n)^0.5)

TK2=transform(TK2,ET=E/2)

# **** 1991 *******
d=subset(TK2,year==1991,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
# **** dT ******
dd=d[!is.na(d$VdT),]
names(dd)[2]="SF"
MM2<-nls(k~ka*(SF)/((S50a)+(SF))+kb*(SF)/((S50b)+(SF)), 
    start=list(ka=.15,kb=.15,S50a=1,S50b=15),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM2))
plotCol(MM2,dd,"dT")

names(dd)[2]="T"
MM2<-nls(k~(kt1*2*T/KT+2*kt2*T^2/(KT*KTT))/
        (1+2*T/KT+T^2/(KT*KTT)), 
    start=list(kt1=.1,kt2=.25,KT=0.5,KTT=15),
#    data=dd,weights=rep(1,length(k)),trace=TRUE)
    data=dd,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(summary(MM2))



# **** dC ******
dd=d[!is.na(d$VdC),]
names(dd)[3]="SF"
MM<-nls(k~ka*(SF)/((S50a)+(SF)),start=list(ka=.1,S50a=3),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM))
plotCol(MM,dd,"dC")


# **** 1999 *******
d=subset(TK2,year==1999,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
# **** dT ******
dd=d[!is.na(d$VdT),]
names(dd)[2]="SF"
MM2<-nls(k~ka*(SF)/((S50a)+(SF))+kb*(SF)/((S50b)+(SF)), 
    start=list(ka=.15,kb=.15,S50a=1,S50b=15),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM2))
plotCol(MM2,dd,"dT")

names(dd)[2]="T"
MM2<-nls(k~(kt1*2*T/KT+2*kt2*T^2/(KT*KTT))/
        (1+2*T/KT+T^2/(KT*KTT)), 
    start=list(kt1=.1,kt2=.25,KT=0.5,KTT=15),
#    data=dd,weights=rep(1,length(k)),trace=TRUE)
    data=dd,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(summary(MM2))


# **** dC ******
dd=d[!is.na(d$VdC),]
names(dd)[3]="SF"
MM<-nls(k~ka*(SF)/((S50a)+(SF)),start=list(ka=.1,S50a=3),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
#MM2<-nls(k~ka*(SF)/((S50a)+(SF))+kb*(SF)/((S50b)+(SF)),start=list(ka=.1,kb=.3,S50a=15,S50b=25),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM))
plotCol(MM,dd,"dC")

# **** 2003 *******
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(frstAut=="Wang")&(seq=="wt"),select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
# **** dT ******
dd=d[!is.na(d$VdT)&(d$dC==0),]
names(dd)[2]="SF"
MM2<-nls(k~ka*(SF)/((S50a)+(SF))+kb*(SF)/((S50b)+(SF)), 
    start=list(ka=.15,kb=.15,S50a=1,S50b=15),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM2))
plotCol(MM2,dd,"dT")
# **** dC ******
dd=d[!is.na(d$VdC)&(d$dT==0),]
names(dd)[3]="SF"
MM<-nls(k~ka*(SF)/((S50a)+(SF)),start=list(ka=.1,S50a=3),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
#MM2<-nls(k~ka*(SF)/((S50a)+(SF))+kb*(SF)/((S50b)+(SF)),start=list(ka=.1,kb=.3,S50a=3,S50b=8),data=dd,weights=rep(1,length(k)),algorithm="port",lower=0)
print(summary(MM))
plotCol(MM,dd,"dC")





