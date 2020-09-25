rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
#TK2=transform(TK2,ET=E/2)
#TK2
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
#d=subset(d,!(((T==1)|(T==2))&(t==0)&(c==1)&(It==1)))  
d=subset(d,!((T>100)&(t==10)&(c==0)&(It==1))) 
#d=subset(d,!((t==0)&(c==0)&(C==0)&(It==1))) # taking out control instead of t=1 takes SSE to .075, instead of .03 so
                                            # it really is t=1 that has to go, and not t=0
d

source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",P=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.A=0,kt.B=0,kc.A="kt.A",kc.B="kt.A",
          KAT=0,KBT=0,KAC=0,KBC="KAC",
          KAt=0, KBt=0, KAc=0,
            KBc=0
#          p=0.5
      ))


(id="DEFG:DE.DE.DEFG")
mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It/2
        +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
    start=list(
        kt.A=.2, 
        kt.B=1,
        kc.A=.2,
        kc.B=.6,
        KAT=.46,
        KBT=106,
        KAC=8.5,
        KBC=8.2,
        KAt=1,
        KBt=10,
        KAc=1,
        KBc=.8),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#c(0.233, 0.978, 0.305, 0.613, 1.112, 17.95, 30.84,  6.200    0.003    6.916   13.851    0.624
s<-summary(mod)$parameters[,"Estimate"]
g$ardP=as.list(s)

cd24=ARD12(g$ardP)

##0.53(0.08)  0.81(0.08)  0.50(0.09)  0.37(0.09)  2.91(0.69)  45.62(15.90)  13.36(3.39) 5.22(0.98)
#(id="DEDG:DE.DE.DEFG")
#kt.A=.53 
#kt.B=.8
#kc.A=.5
#kc.B=.4
#KAT=2.9
#KBT=45.6
#KAC=13.4
#KBC=5.22
#KAt=1
#KBt=1
#KAc=13
#KBc=.6
#mod<-nls(k~((.53*T/2.9/(1+T/2.9+C/13.4+t/KAt+c/KAc))
#          +(.8*T/45.6/(1+T/45.6+C/5.2+t/KBt+c/KBc)))*It/2
#        +((.5*C/13.4/(1+T/2.9+C/13.4+t/KAt+c/KAc))
#          +(.4*C/5.2/(1+T/45.6+C/5.2+t/KBt+c/KBc)))*Ic/2, 
#    start=list(#KAt=.3,
#KBt=1,
#KAc=1,
#KBc=.99
#),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmodAll(mod,d,tab,id)
#s<-summary(mod)$parameters[,"Estimate"]
##s[1:length(s)]=rep(1,length(s))
#ardP=as.list(s)
#ardP["kc.A"]=ardP["kt.A"]
#ardP
#ard2dim(ardP)
#
#(id="DEXG:DE.DE.DEIG")
#mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+t/KAt+c/KAc))
#          +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It/2
#        +(
#          (kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
#    start=list(
#        kt.A=.2, 
#        kt.B=1,
#        #        kc.A=.2,
#        kc.B=.6,
#        KAT=.46,
#        KBT=106,
##        KAC=8.5,
#        KBC=8.2,
#        KAt=1,
#        KBt=10,
#        KAc=1,
#        KBc=.8),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#s<-summary(mod)$parameters[,"Estimate"]
##s[1:length(s)]=rep(1,length(s))
#ardP=as.list(s)
#ardP["kc.A"]=ardP["kt.A"]
#ardP["KAC"]=Inf
#ardP
#ARD12(ardP)
#
#
#s<-summary(mod)$parameters[,"Estimate"]
#save(s,file="C:/Users/radivot/case/active/papers/TK2/ARDall.RData")

library(hwriter)
hwrite(t(tab),"C:/Users/radivot/case/active/papers/TK2/ARDall.html",cellspacing=0,cellpadding=2)

#KAt=0.02 # (Barosso 2003 referenced to 1992, tight value)
#(id="DEFG:DEfGKLMN")
#mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It/2
#        +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
#    start=list(
#        kt.A=.2, 
#        kc.A=.2,
#        kt.B=1,
#        kc.B=.6,
#        KAT=.46,
#        KAC=8.5,
##        KAt=1,
#        KAc=1,
#        KBT=106,
#        KBC=8.2,
#        KBt=10,
#        KBc=.8),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)



#(id="DEFG:DEFGKEMN")
#mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kt.B*T/KBT/(1+T/KBT+C/KAC+t/KBt+c/KBc)))*It/2
#        +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kc.B*C/KAC/(1+T/KBT+C/KAC+t/KBt+c/KBc)))*Ic/2, 
#    start=list(
#        kt.A=.2, 
#        kc.A=.2,
#        kt.B=1,
#        kc.B=.6,
#        KAT=.46,
#        KAC=8.5,
#        KAt=1,
#        KAc=1,
#        KBT=106,
##        KBC=8.2,
#        KBt=10,
#        KBc=.8),
##    start=list(
##        kt.A=.2, kt.B=1,
#    ##        kc.A=.2,
##        kc.B=.6,
##        KAT=.46,KBT=106,
##        KAC=8.5, #KBC=8.2,
##        KAt=1,KBt=10,
##        KAc=1,KBc=.8),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)

