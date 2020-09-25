rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
d=subset(TK2,(year>=2003)&(seq=="wt"))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[1:2]=c("T","C")
d1=subset(d,(((It==1)&(C==0))|((Ic==1)&(T==0)))&(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(frstAut=="Wang"))
d2=subset(d,(year==2003)&(frstAut=="Barroso")&(state=="dimer"))
d3=subset(d,(year==2003)&(frstAut=="Barroso")&(state=="tetramer"))
d4=subset(d,(year==2005)&(state=="dimer"))
d5=subset(d,(year==2005)&(state=="tetramer"))
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",P=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",
          kc.C="kt.T",kc.CC="X",
          KT=0,KTT=0,KC=0,KCC=Inf))

Xg = c(0.05,0.075, 0.1, 0.2,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
mods=NULL
df=d1
i=1
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.9,KT=13,KTT=160,KC=18),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDX:DE.DI")
mod$data="2003 Wang"
mod$d=df
mods[[i]]=mod
tab=plotmod(mod,tab,kill=T)

df=d2
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.8,KT=32,KTT=67,KC=35),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDX:DE.DI")
mod$data="2003 Barroso Dimer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

df=d3
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.2,KT=32,KTT=67,KC=35),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDX:DE.DI")
mod$data="2003 Barroso Tetramer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

df=d4
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.8,kc.C=.3,KT=.3,KTT=5,KC=5),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDFX:DE.DI")
mod$data="2005 Barroso Dimer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)


df=d4
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC +2*kc.C*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.3,kt.TT=.53,kc.C=.35,KT=.2,KTT=4.6,KC=4.5),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFF:DE.DD")
mod$data="2005 Barroso Dimer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.4,KT=.3,KTT=5,KC=5),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDFX:DE.DI")
mod$data="2005 Barroso Dimer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

df=d5
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC +2*kc.C*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.3,kt.TT=.53,kc.C=.35,KT=.2,KTT=4.6,KC=4.5),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFF:DE.DD")
mod$data="2005 Barroso Tetramer"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)
tab

tabs=mkTab(mods)
names(tabs)[2]<-"Model"
names(tabs)[1]<-"Data"
tabs[,"kc.CC"]="X"
tabs[c(5,7),"kc.CC"]="kc.C"
tabs[,"KCC"]="Inf"
tabs[c(5,7),"KCC"]="KC"
library(hwriter)
hwrite(tabs,"C:/Users/radivot/case/active/papers/TK2/barroso.html",cellspacing=0,cellpadding=2)

aics=as.numeric(tabs$AIC)
mins=tapply(aics,tabs$Data,min)
picks=NULL
for (i in 1:length(mins)) picks[i]=which(aics==mins[i])
picks=sort(picks)
picks
plotmods(mods,d,picks,kill=T,Xg=Xg)

