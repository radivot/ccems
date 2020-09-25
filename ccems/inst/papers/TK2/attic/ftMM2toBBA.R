library(ccems)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 12, height = 8,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfrow=c(3,4),mar=c(4,4,2,1)+.1)
for (j in 1:11){
  d=subset(TK2BBA77,donor==j)
  print(d)
  MMmod<-nls(V~ka*S/(S50a+S)+kb*S/(S50b+S), 
      start=list(ka=10,kb=200,S50a=1,S50b=50),data=d,weights=V/V,algorithm="port",lower=0)
#  start=list(ka=.15,kb=.15,S50a=1,S50b=15),data=d,weights=V/V,algorithm="port",lower=0)
  print(summary(MMmod))
  plot(d$S,d$V,xlab="Total [dT]",log="xy", ylab="k (1/sec) ", 
      main=paste("Donor #",d[1,"donor"],sep="") )
  lgx=log(d$S)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  lines(fineX,predict(MMmod,list(S=fineX)),col="black",lwd=1)
  mtext(paste("S50a = ",format(MMmod$m$getPars()["S50a"],digits=3),sep=""),line=-2.2,side=1,font=1,cex=0.7)
  mtext(paste("S50b = ",format(MMmod$m$getPars()["S50b"],digits=3),sep=""),line=-1.2,side=1,font=1,cex=0.7)
}
