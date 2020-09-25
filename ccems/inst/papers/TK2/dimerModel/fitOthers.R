rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
DNA="H121N"
DNA="wt"
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[1:2]=c("T","C")
dtm=subset(d,(It==1)&(C==0))
dcm=subset(d,(Ic==1)&(T==0))
dm=rbind(dtm,dcm) # m for marginals
d1=subset(dm,year==1991)
d2=subset(dm,year==1999)
d3=subset(dm,(year==2003)&(frstAut=="Wang"))
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",P=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",
          kc.C="kt.T",kc.CC="X",
          KT=0,KTT=0,KC=0,KCC=Inf))
Xg = c(0.1,0.15,0.2,0.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
mods=NULL



i=1
df=d3
mod<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=10,KTT=100,KC=10),
        data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFF:DE.DD")
mod$data="2003 Wang"
mod$d=df
mods[[i]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.8,kc.C=.6,KT=10,KTT=100,KC=10),
        data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDFF:DE.DD")
mod$data="2003 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.52,KT=8,KTT=7,KC=16),
        data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDD:DE.DD")
mod$data="2003 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC)*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
        start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=10,KTT=100,KC=10),
        data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFX:DE.DI")
mod$data="2003 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
        start=list(kt.T=.4,KT=10,KTT=100,KC=10),
        data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDX:DE.DI")
mod$data="2003 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)


df=d1
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=.3,KTT=16,KC=5),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFF:DE.DD")
mod$data="1991 Munch-Petersen"
mod$d=df
mods[[i<-i+1]]=mod

mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.2,kc.C=.2,KT=.3,KTT=23,KC=43),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDFF:DE.DD")
mod$data="1991 Munch-Petersen"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,KT=.6,KTT=23,KC=43),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDD:DE.DD")
mod$data="1991 Munch-Petersen"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)



df=d2
mod<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.2,kt.TT=.35,kc.C=.4,KT=0.86,KTT=8.5,KC=20),
    data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEFF:DE.DD")
mod$data="1999 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kt.TT*2*C/KC+2*kt.TT*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.2,kt.TT=.35,KT=0.86,KTT=8.5,KC=20),
        data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DEEE:DE.DD")
mod$data="1999 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

mod<-
    nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.35,KT=0.86,KTT=8.5,KC=20),
        data=df,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
(mod$id="DDDD:DE.DD")
mod$data="1999 Wang"
mod$d=df
mods[[i<-i+1]]=mod
tab=plotmod(mod,tab,Xg=Xg)

tabs=mkTab(mods)
names(tabs)[2]<-"Model"
names(tabs)[1]<-"Data"
tabs[,"kc.CC"]="kc.C"
tabs[c(4,5),"kc.CC"]="X"
tabs[c(3,8,11),"kc.CC"]="kt.T"
tabs[10,c("kc.C","kc.CC")]="kt.TT"
tabs[,"KCC"]="KC"
tabs[c(4,5),"KCC"]="Inf"
tabs

library(hwriter)
hwrite(tabs,"C:/Users/radivot/case/active/papers/TK2/others.html",cellspacing=0,cellpadding=2)

aics=as.numeric(tabs$AIC)
mins=tapply(aics,tabs$Data,min)
picks=NULL
for (i in 1:length(mins)) picks[i]=which(aics==mins[i])
picks=sort(picks)
picks
plotmods(mods,d,picks,kill=T,Xg=Xg)

