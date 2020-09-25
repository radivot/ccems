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
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(frstAut=="Wang")&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 

showBad=TRUE
showBad=FALSE
if (!showBad) d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
#d=subset(d,!((C==5)&(T==161))) # kill outlier with higher rate than no inhibitor 
d
source("common.R")

options(stringsAsFactors = FALSE)
############# START WITH 6 IN TABLE BY LOGICAL PROGRESSION ###############
(tab=data.frame(id="XX",P=rep(0,14),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",kt.TC=c(rep("kt.T",14)),
          kc.C="kt.T",kc.CC=c(rep("kc.C",6),rep("X",8)),kc.TC=c(rep("kc.C",6),rep("kt.TC",8)),
          KT=0,KTT=0,KC=0,KCC=c(rep("KC",6),rep(Inf,8)),KTC="KC"))

(id="DEFKKM:DE.DEF")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KCC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KCC))/2, 
start=list(kt.T=.66,kt.TT=.83,kt.TC=0.63,
           kc.C=.43,kc.TC=0.33,
           KT=6,KTT=56,
           KC=9,KCC=8,KTC=9),
#    start=list(kt.T=.8,kt.TT=.9,kt.TC=0.25,kc.C=.8,kc.TC=.61,KT=8.5,KTT=90,KC=19,KCC=800,KTC=22), # unstable
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=T)

(id="DEFKKM:DE.DDF")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,
               kc.C=.6,kc.TC=0.2,
               KT=5,KTT=50,
               KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEFKKM:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,
               kc.C=.6,kc.TC=0.2,
               KT=5,KTT=50,
               KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEDKKM:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,
               kc.C=.6,kc.TC=0.2,
               KT=5,KTT=50,
               KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEDKKK:DE.DDD")
KTic=ifelse(DNA=="wt",5,.1)
KTTic=ifelse(DNA=="wt",50,7)
kt.Tic=ifelse(DNA=="wt",.4,.1)
kc.Cic=ifelse(DNA=="wt",.6,.06)
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+
          (kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=kt.Tic,kt.TT=.8,
               kc.C=kc.Cic,
               KT=KTic,KTT=KTTic,
               KC=6),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDDKKK:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=.4,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEFKXM:DE.DIF")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It
          +(kc.C*2*C/KC+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC)/2, 
    control=list(maxiter=500),
    start=list(
        kt.T=.9,kt.TT=.9,kt.TC=0.63,
        kc.C=.83,kc.TC=0.6,
        KT=10,KTT=100,
        KC=18,KTC=20),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDFKXM:DE.DIF")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It
          +(kc.C*2*C/KC+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC)/2, 
    control=list(maxiter=500),
    start=list(
        kt.T=.9,kt.TC=0.63,
        kc.C=.83,kc.TC=0.6,
        KT=10,KTT=100,
        KC=18,KTC=20),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDFKXM:DE.DID")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It
          +(kc.C*2*C/KC+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC)/2, 
    start=list(
        kt.T=.66,kt.TC=0.63,
        kc.C=0.8,kc.TC=0.33,
        KT=6,KTT=56,
        KC=18),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDFDXM:DE.DIF")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It
          +(kt.T*2*C/KC+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC)/2, 
    start=list(
        kt.T=.8,kt.TC=0.63,
        kc.TC=0.33,
        KT=6,KTT=56,
        KC=18,KTC=19),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDFDXM:DE.DID")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It
          +(kt.T*2*C/KC+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC)/2, 
    start=list(
        kt.T=.66,kt.TC=0.63,
        kc.TC=0.33,
        KT=6,KTT=56,
        KC=18),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)


(id="DDFDXF:DE.DID")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It
          +(kt.T*2*C/KC+kt.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC)/2, 
    start=list(
        kt.T=.66,kt.TC=0.53,
#        kc.TC=0.53,
        KT=6,KTT=56,
        KC=18),
#    data=d,weights=1/k^2,trace=TRUE) # resids at low [dT] still there
data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DFFDXF:DE.DID")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TC*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It
          +(kt.T*2*C/KC+kt.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC)/2, 
    control=list(maxiter=500),
    start=list(
        kt.T=.66,kt.TC=0.53,
#        kc.TC=0.53,
        KT=6,KTT=56,
        KC=18),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)
tab[13,"kt.TT"]="kt.TC"

(id="DDDDXM:DE.DID")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It
          +(kt.T*2*C/KC+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC)/2, 
    start=list(
        kt.T=.66,#kt.TC=0.53,
        kc.TC=0.53,
        KT=6,KTT=56,
        KC=18),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDXDXX:DE.DII")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.66,
#        kt.TT=.83,
#        kc.C=.43,
        KT=6,KTT=56,
        KC=9),
#    start=list(kt.T=.8,kt.TT=.9,kt.TC=0.25,kc.C=.8,kc.TC=.61,KT=8.5,KTT=90,KC=19,KCC=800,KTC=22), # unstable
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)


tab1=tab
tab1
names(tab1)[1]<-"Model"
#names(tab1)[11:15]<-c("KD_T","KDT_T","KD_C","KDC_C","KDT_C")
library(hwriter)
hwrite(tab1,"C:/Users/radivot/case/active/papers/TK2/tab1.html",cellspacing=0,cellpadding=2)

