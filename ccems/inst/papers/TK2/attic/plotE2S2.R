library(ccems)
rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
print(TK2)

## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 3, height = 6,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
par(mfcol=c(2,1),mar=c(4,4,2,1)+.1)
d=subset(TK2,(dT>0)&(year==1999),select=c(dT,E,k,frstAut,year))
names(d)[1:2]=c("ST","ET")
print(d)
topology <- list(  
    heads=c("E1S0","E2S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site 
            m=c("E1S1"),
            d= c("E2S1","E2S2")  
        ) # d for dimer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=T)
#g
plot(d$ST,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec) ", 
      main=paste(d[1,"frstAut"],d[1,"year"]))

  
  lgx=log(d$ST)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  lines(fineX,predict(MMmod,list(TT=fineX)),col="black",lwd=1)
  plot(MMmod$m$fitted(),MMmod$m$resid(),xlab="Fitted Value",
      ylab="Relative Residual",mar=c(2,2,0,1)+.1)


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



