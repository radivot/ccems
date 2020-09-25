rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,ET=E/2)
#TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(seq==DNA),select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[2:3]=c("T","C")
d=transform(d,It=as.numeric(!is.na(d$VdT)),Ic=as.numeric(!is.na(d$VdC)))
d=subset(d,!((C==5))) # kill outlier with higher rate than no inhibitor 
#d=subset(d,!((C==5)&(T==161))) # kill outlier with higher rate than no inhibitor 
d

aics=NULL

plotmod<-function(mod,d,tab,id,jj=1,kill=FALSE,pause=FALSE){
  if (class(mod)=="nls"){
    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
    if (DNA=="wt") {Cg=c(0,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
    pd=rbind(data.frame(T=rep(Xg,length(Cg)),C=rep(Cg,each=length(Xg)),It=1,Ic=0),
        data.frame(T=rep(Tg,each=length(Xg)),C=rep(Xg,length(Tg)),It=0,Ic=1) )
    pd=data.frame(pd,Ek=predict(mod,pd))
    if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
    if (.Platform$OS.type=="windows") 
      windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
    par(mfcol=c(1,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,3,1))
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
    for (i in 1:length(Cg)) 
    {
      with(d[(d$It>0)&(d$C==Cg[i]),],points(T,k,col=i,pch=i))
      with(pd[(pd$It>0)&(pd$C==Cg[i]),],lines(T,Ek,col=i))
    }
    legend("topleft",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
    
    print(s<-summary(mod))
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
    for (i in 1:length(Tg)) 
    {
      with(d[(d$Ic>0)&(d$T==Tg[i]),],points(C,k,col=i,pch=i))
      with(pd[(pd$Ic>0)&(pd$T==Tg[i]),],lines(C,Ek,col=i))
    }
    legend("topleft",legend=paste("[dT]=",Tg,sep=""),pch=1:length(Tg),col=1:length(Tg),bty="n")
    (p=dim(s$parameters)[1])
    ed=data.frame(d,Ek=predict(mod,d))
    SSE=sum((ed$Ek-ed$k)^2)
    N=dim(d)[1]
    P=p+1
    sig=format(sqrt(SSE/(N-p)),digits=2)
    AICc=format(N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N,digits=4)
    SSE=format(SSE,digits=3)
#    AICuc=format(N*log(SSE/N)+2*P + N*log(2*pi) + N,digits=4)
#    AICchk=format(AIC(mod),digits=4)
    title(paste(id," (",DNA,", AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=1.7,cex.main=1)
    title(paste(sub(" ","",names(s$parameters[,1])),"=",format(s$parameters[,1],digits=1,trim=T),sep="",collapse=";  "),
        outer=TRUE,line=0.7,cex.main=0.8)
    if (pause) locator(n=1)
    par(mfrow=c(1,1))
    vals=paste(format(s$parameters[,1],digits=2,trim=T),"(",format(s$parameters[,2],digits=1,trim=T),")",sep="")
    names(vals)<-names(s$parameters[,1])
    tmp=c(id=id,p=p,SSE=SSE,AIC=AICc,vals)
#tmp=c(id=id,p=p,aic=format(AIC(mod),digits=4),format(s$parameters[,1],digits=2))
    if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
    tab[jj,names(tmp)]=tmp
    } else {tab[jj,1]=id; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  tab
}

options(stringsAsFactors = FALSE)
############# START WITH 6 IN TABLE BY LOGICAL PROGRESSION ###############
(tab=data.frame(id="XX",p=rep(0,13),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",kt.TC=c(rep("kt.T",13)),
          kc.C="kt.T",kc.CC=c(rep("kc.C",6),rep("X",7)),kc.TC=c(rep("kc.C",6),rep("kt.TC",7)),
          KT=0,KTT=0,KC=0,KCC=c(rep("KC",6),rep(Inf,7)),KTC="KC"))

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
tab[12,"kt.TT"]="kt.TC"

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

tab1=tab
tab1
names(tab1)[1]<-"Model"
library(hwriter)
hwrite(tab1,"C:/Users/radivot/case/active/papers/TK2/tab1.html",cellspacing=0,cellpadding=2)

############# FINISH 6 IN TABLE BY LOGICAL PROGRESSION ###############
#################################  new tabe with 24 entries, brute force #########
(tab=data.frame(id="XX",p=rep(0,24),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",kt.TC="kt.T",kc.C=0,kc.CC="kc.C",kc.TC="kc.C",KT=0,KTT=0,KC=0,KCC="KC",KTC="KC"))
(id="DEF.DDE:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,kc.TC=0.2,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
(id="DEF.DDE:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,kc.TC=0.2,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DEF.DED:DE.DDE")
e=try(mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,kc.CC=0.2,KT=5,KTT=40,KC=9,KTC=7),
    control=list(maxiter=2000),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE))
if (class(e)=="try-error") 
    tab=plotmod("failed",d,tab,id,jj=3) 
(id="DEF.DED:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,kc.CC=0.2,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
tab
####
(id="DEF.DDD:DE.DDE")
e=try(mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,KT=5,KTT=50,KC=9,KTC=11),
        control=list(maxiter=2000),
        data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE))
if (class(e)=="try-error") 
  tab=plotmod("failed",d,tab,id,jj=3) else tab=plotmod(mod,d,tab,id) 
(id="DEF.DDD:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kt.TC=0.2,kc.C=.6,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
tab
######################################################
(id="DED.DDE:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,kc.TC=0.2,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
(id="DED.DDE:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,kc.TC=0.2,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DED.DED:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=0.2,kc.CC=0.2,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id,7)
(id="DED.DED:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,kc.CC=0.2,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DED.DDD:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=0.2,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id,7)
(id="DED.DDD:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
tab
######################################################
(id="DDE.DDE:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC)), 
    start=list(kt.T=.6,kt.TC=.4,kc.C=.4,kc.TC=0.3,KT=5,KTT=50,KC=9,KTC=11)/2,
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
(id="DDE.DDE:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
    start=list(kt.T=.6,kt.TC=.4,kc.C=.4,kc.TC=0.3,KT=5,KTT=50,KC=9)/2,
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DDE.DED:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.6,kt.TC=.4,kc.C=0.4,kc.CC=0.3,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
(id="DDE.DED:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.6,kt.TC=.4,kc.C=.4,kc.CC=0.3,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DDE.DDD:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.6,kt.TC=.4,kc.C=0.4,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
(id="DDE.DDD:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.6,kt.TC=.4,kc.C=.4,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
tab

######################################################
(id="DDD.DDE:DE.DDE")
e=try(mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.7,kc.C=.4,kc.TC=0.3,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE))
if (class(e)=="try-error") 
  tab=plotmod("failed",d,tab,id,jj=13) else tab=plotmod(mod,d,tab,id) 
#tab=plotmod(mod,d,tab,id)
(id="DDD.DDE:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=.4,kc.TC=0.3,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
####
(id="DDD.DED:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=0.4,kc.CC=0.3,KT=5,KTT=50,KC=9,KTC=11),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id,7)
(id="DDD.DED:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=.4,kc.CC=0.3,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
print(tab)
####
(id="DDD.DDD:DE.DDE")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=0.4,KT=5,KTT=40,KC=9,KTC=7),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id,7)
(id="DDD.DDD:DE.DDD")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=.4,KT=5,KTT=50,KC=9),
    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
tab=plotmod(mod,d,tab,id)
print(tab)
hwrite(tab,"C:/Users/radivot/case/active/papers/TK2/tab2.html")

######################################################
#######################################################
#(id="DXE.DDE:DI.DDE")
#mod<-nls(k~((kt.T*2*T/KT+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.3,kt.TC=.2,kc.C=.2,kc.TC=0.15,KT=5,KC=9,KTC=11),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#(id="DXE.DDE:DI.DDD")
#mod<-nls(k~((kt.T*2*T/KT+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.3,kt.TC=.2,kc.C=.2,kc.TC=0.15,KT=5,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#####
#(id="DXE.DED:DI.DDE")
#mod<-nls(k~((kt.T*2*T/KT+kt.TC*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.3,kt.TC=.2,kc.C=0.2,kc.CC=0.15,KT=5,KC=9,KTC=11),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id,7)
#(id="DXE.DED:DI.DDD")
#mod<-nls(k~((kt.T*2*T/KT+kt.TC*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.3,kt.TC=.2,kc.C=.2,kc.CC=0.15,KT=5,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#tab
#(id="DXD.DDE:DI.DDE")
#e=try(mod<-nls(k~((kt.T*2*T/KT+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KTC))*Ic)/
#            (1+2*T/KT+2*T*C/(KT*KTC)+2*C/KC+2*C/KC+C^2/(KC*KC)), 
#        start=list(kt.T=.35,kc.C=.2,kc.TC=0.15,KT=5,KC=9,KTC=11),
#        data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE))
#if (class(e)=="try-error") 
#  tab=plotmod("failed",d,tab,id,jj=21) else tab=plotmod(mod,d,tab,id) 
##tab=plotmod(mod,d,tab,id)
#(id="DXD.DDE:DI.DDD")
#mod<-nls(k~((kt.T*2*T/KT+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.TC*2*T*C/(KT*KC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.4,kc.C=.2,kc.TC=0.15,KT=5,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#####
#(id="DXD.DED:DI.DDE")
#mod<-nls(k~((kt.T*2*T/KT+kt.T*2*T*C/(KT*KTC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KTC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KTC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.4,kc.C=0.2,kc.CC=0.15,KT=5,KC=9,KTC=11),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id,7)
#(id="DXD.DED:DI.DDD")
#mod<-nls(k~((kt.T*2*T/KT+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.CC*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
#        (1+2*T/KT+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.4,kc.C=.2,kc.CC=0.15,KT=5,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)


#(id="DED.DDD:DE.DDD")
#mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC))*Ic)/
#        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.2,kt.TT=.4,kc.C=.3,KT=5,KTT=50,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#

#
## DDX.DDX.DDI.DD
#mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
#        (1+2*T/KT+T^2/(KT*KT)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.2,kc.C=.3,KT=5,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#AIC(mod)
#plotmod(mod,d,kill=FALSE)


## ET_C=Inf
## DEX.DDX.DEI.DD
#mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
#        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.2,kt.TT=.2,kc.C=.3,KT=5,KTT=8,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
## convergence failure on kt.TT/KTT
#
## DDX.DDX.DEI.DD
#mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
#        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC)), 
#    start=list(kt.T=.2,kc.C=.3,KT=5,KTT=8,KC=9),
#    data=d,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#print(summary(mod))
#AIC(mod)
#plotmod(mod,d)
## bad fit => skip the KTC=Inf hypothesis

################


#with(pm[pm$It>0,],lines(T,Ek,col="violet"))
#with(pm[pm$Ic>0,],lines(C,Ek,col="black"))
#legend("topleft",legend=c("dT kinase activity","dC kinase activity"),pch=1,col=c("violet","black"),bty="n")
#ph=pd[(pd$T!=0)&(pd$C!=0),]
#ph
#dh=d[(d$T!=0)&(d$C!=0),]
#dh
#
#
#nT=sum(dh$It)
#nC=sum(dh$Ic)
#Xgh <- c(0.5,1,2,5,10,20,50,100,200)
#library(rgl)
##with(dh,plot3d(log10(T),log10(C),k,col=c(rep("violet",nT),rep("blue ",nC)),
##type="s",radius=c(rep(.05,nT+nC)),ylab="[dC]",xlab="[dT]",zlab="1/sec"))
#with(dh[dh$It==1,],plot3d(log10(T),log10(C),k,col="violet",
#type="s",radius=.05,ylab="[dC]",xlab="[dT]",zlab="1/sec"))
##material3d(alpha=1)
#surface3d(log10(Xgh),log10(Xgh),alpha=0.7,t(matrix(ph$Ek[ph$It==1],ncol=length(Xgh))),col="violet")
#
#with(dh[dh$Ic==1,],plot3d(log10(T),log10(C),k,col="violet",
#type="s",radius=.05,ylab="[dC]",xlab="[dT]",zlab="1/sec"))
##material3d(alpha=1)
#surface3d(log10(Xgh),log10(Xgh),alpha=0.7,t(matrix(ph$Ek[ph$Ic==1],ncol=length(Xgh))),col="violet")
## the snipping tool in vista was used to capture the images in Figure 1

