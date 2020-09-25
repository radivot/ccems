rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))


DNA="wt"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
dtw=subset(d,(It==1)&(C==0)&(c==0)&(t==0))
dtw
Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 4, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(1,1),mar=c(4.2,4,.5,1)+.1,oma=c(0,0,0,0))
plot(Xg,type="n",ylim=range(c(0,dtw$k)),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
pd=data.frame(T=Xg,It=1)

mod<-nls(k~(kt1*2*T/KT+2*kt1*T^2/(KT*KTT))/(1+2*T/KT+T^2/(KT*KTT)), 
    start=list(kt1=.4,KT=5,KTT=50),
    data=dtw,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(s<-summary(mod))
pdw=data.frame(pd,Ek=predict(mod,pd))
with(dtw,points(T,k))
with(pdw,lines(T,Ek))


DNA="H121N"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
dtm=subset(d,(It==1)&(C==0)&(c==0)&(t==0))
dtm
KT=10;KTT=75
mod<-nls(k~(kt1*2*T/KT+2*kt1*T^2/(KT*KTT))/(1+2*T/KT+T^2/(KT*KTT)), 
    start=list(kt1=.2),
    data=dtm,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(s<-summary(mod))
pdm=data.frame(pd,Ek=predict(mod,pd))
with(dtm,points(T,k,pch=2))
with(pdm,lines(T,Ek))
ed=data.frame(dtw,Ek=predict(mod,dtw))
SSE=sum((ed$Ek-ed$k)^2)
SSE


mod<-nls(k~kt1*(T/KT)^n/(1+(T/KT)^n), 
    start=list(kt1=.2,KT=10,n=1),
    data=dtm,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(s<-summary(mod))
#pdm=data.frame(pd,Ek=predict(mod,pd))
#      with(pdm,lines(T,Ek,col="red"))

mod<-nls(k~kt1*(T/KT)/(1+(T/KT)), 
    start=list(kt1=.2,KT=10),
    data=dtm,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(s<-summary(mod))
pdm=data.frame(pd,Ek=predict(mod,pd))
with(pdm,lines(T,Ek,col="red"))
legend("topleft",legend=c("Wild Type","H121N"),pch=c(1,2),bty="n")
ed=data.frame(dtw,Ek=predict(mod,dtw))
SSE=sum((ed$Ek-ed$k)^2)
SSE



KTT=Inf;kt1=0.4
mod<-nls(k~(kt1*2*T/KT+2*kt1*T^2/(KT*KTT))/(1+2*T/KT+T^2/(KT*KTT)), 
    start=list(KT=10),
    data=dtm,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
print(s<-summary(mod))
pdm=data.frame(pd,Ek=predict(mod,pd))
with(pdm,lines(T,Ek,col="blue"))
#
#KT=10;kt1=0.4
#mod<-nls(k~(kt1*2*T/KT+2*kt1*T^2/(KT*KTT))/(1+2*T/KT+T^2/(KT*KTT)), 
#    start=list(KTT=50),
#    data=dtm,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)



