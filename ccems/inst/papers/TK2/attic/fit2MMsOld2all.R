library(ccems)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 12, height = 6,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(2,4),mar=c(4,4,2,1)+.1)
for (j in c(1991,1999,2003)){
  d=subset(TK2,(fg==ifelse(j==1999,2,4))&(dCTP==0)&(dTTP==0)&(dC==0)&(dT>0)&(year==j),select=c(dT,E,k,frstAut,year))
  names(d)[1:2]=c("TT","ET")
  print(d)
  MMmod<-nls(k~ka*(TT)/((S50a)+(TT))+kb*(TT)/((S50b)+(TT)), 
      start=list(ka=.15,kb=.15,S50a=1,S50b=15),data=d,weights=1/k^2,algorithm="port",lower=0)
  print(summary(MMmod))
  plot(d$TT,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec) ", 
      main=paste(d[1,"frstAut"],d[1,"year"]))
  lgx=log(d$TT)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  lines(fineX,predict(MMmod,list(TT=fineX)),col="black",lwd=1)
  plot(MMmod$m$fitted(),MMmod$m$resid(),xlab="Fitted Value",
      ylab="Relative Residual",mar=c(2,2,0,1)+.1)
}

j=1999
d=subset(TK2,(fg==ifelse(j==1999,2,4))&(dCTP==0)&(dTTP==0)&(dC==0)&(dT>0)&(dT<1.2)&(year==j),select=c(dT,E,k,frstAut,year))
names(d)[1:2]=c("TT","ET")
MMmod<-nls(k~ka*(TT)/(S50a+TT), 
    start=list(ka=.15,S50a=1),data=d,weights=1/k^2,algorithm="port",lower=0)
print(summary(MMmod))
plot(d$TT,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec) ", 
    main=paste(d[1,"frstAut"],d[1,"year"]))
lgx=log(d$TT)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
lines(fineX,predict(MMmod,list(TT=fineX)),col="black",lwd=1)
plot(MMmod$m$fitted(),MMmod$m$resid(),xlab="Fitted Value",
    ylab="Relative Residual",mar=c(2,2,0,1)+.1)



