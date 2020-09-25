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
#TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
#d=subset(d,!(((T==1)|(T==2))&(t==0)&(c==1)&(It==1)))  
#d

d=subset(d,!((T>100)&(t==10)&(c==0)&(It==1))) 

#d=subset(d,!((t==0)&(c==0)&(C==0)&(It==1))) # taking out control instead of t=1 takes SSE to .075, instead of .03 so
# it really is t=1 that has to go, and not t=0
#d
g=NULL
g$d=d
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",P=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.A=0,kt.B=0,kc.A="kt.A",kc.B="kt.A",
          KAT=0,KBT=0,KAC=0,KBC="KAC",
          KAt=0, KBt=0, KAc=0, KBc=0
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
ardP=as.list(s)
cd24=ARD12(ardP)
library(hwriter)
cd24df=as.data.frame(lapply(cd24,format,digits=2))
hwrite(cd24df,"C:/Users/radivot/case/active/papers/TK2/ARDall.html",cellspacing=0,cellpadding=2)
#cd24df

g$CD=FALSE
g$CD=TRUE
if (g$CD) {  g$params=data.frame(initial=unlist(cd24),final=unlist(cd24),opt=TRUE,constr="none")
}  else {  g$params=data.frame(initial=s,final=s,opt=TRUE,constr="none") }# else refit the ARD model
gg=getEk(g,init=TRUE)

mods=NULL
mod=NULL
for (i in 1:24) {
  g=gg
  g$params[i,"opt"]=FALSE
  e=try(g<-fitModel(g))
  if (class(e)=="try-error") mod="failed" else {
    mod$params=g$params
    mod$SSE=g$SSE
    mod$AIC=g$AIC
  }
  mods[[i]]=mod
}

mods


options(stringsAsFactors = FALSE)
############# START WITH 6 IN TABLE BY LOGICAL PROGRESSION ###############
(tab=data.frame(id="XX",P=0,SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kt.TT="kt.T",kt.TC=0,
          kc.C="kt.T",kc.CC=rep("X",5),kc.TC=rep("kt.TC",5),
          KT=0,KTT=0,KC=0,KCC=Inf,KTC="KC",
          kt.tT="kt.TC",kt.cT=c("kt.tT",rep("kt.TC",4)),
          kc.cC="X",kc.tC=c("kt.tT",rep("kt.TC",4)),
          Kt=0,Kc="Kt",
          KtT=c(0,0,0,0,"KC"),
          KcT=c(0,0,0,"KtT","KC"),
          KtC=c(0,0,0,0,0),
          KcC=Inf,
          Ktt=Inf,
          Ktc=Inf,
          Kcc=Inf))
#kt.T kt.TT kt.TC kt.tT kt.cT kc.C kc.CC kc.TC kc.tC kc.cC KT   KTT KC  KCC KTC KtT KcT KtC  KcC  Kt  Kc   Ktt Ktc Kcc
#0.28   0.6  0.24  0.97  0.23 0.56  0.46  0.61  0.61  0.36 2.1 9.5 10  19  6.5  18 1.2  6.2  26 0.006 1.2 3.5 0.62 7.2

(id="DEFGK,LMNRS:DE.DEF.DEFG.DEFGK")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.cT*2*c*T/(Kc*KcT)
            )*It
          +((kc.C*2*C/KC)+2*kc.CC*C^2/(KC*KCC)+kc.TC*2*T*C/(KT*KTC) 
            + kc.tC*2*t*C/(Kt*KtC)
            + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
        + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
          )/2, 
    control=list(minFactor=1e-16),
    start=cd24,
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)

(id="DEFGK,LMEES:DE.DEF.DEFG.DEFGK")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.cT*2*c*T/(Kc*KcT)
            )*It
          +((kc.C*2*C/KC)+2*kc.CC*C^2/(KC*KCC)+kt.TT*2*T*C/(KT*KTC) 
            + kt.TT*2*t*C/(Kt*KtC)
            + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
          + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
          )/2, 
    control=list(minFactor=1e-16),
    start=cd24[-c(8,9)],
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#kt.T kt.TT kt.TC kt.tT kt.cT kc.C kc.CC kc.TC kc.tC kc.cC KT   KTT KC  KCC KTC KtT KcT KtC  KcC  Kt  Kc   Ktt Ktc Kcc
#0.28   0.6  0.24  0.97  0.23 0.56  0.46  0.61  0.61  0.36 2.1 9.5 10  19  6.5  18 1.2  6.2  26 0.006 1.2 3.5 0.62 7.2

(id="DEDGD,LMEES:DE.DEF.DEFG.DEFGK")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.T*2*c*T/(Kc*KcT)
            )*It
          +((kc.C*2*C/KC)+2*kc.CC*C^2/(KC*KCC)+kt.TT*2*T*C/(KT*KTC) 
            + kt.TT*2*t*C/(Kt*KtC)
            + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
          + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
          )/2, 
    control=list(minFactor=1e-16),
    start=cd24[-c(3,5,8,9)],
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#kt.T kt.TT kt.TC kt.tT kt.cT kc.C kc.CC kc.TC kc.tC kc.cC KT   KTT KC  KCC KTC KtT KcT KtC  KcC  Kt  Kc   Ktt Ktc Kcc
#0.28   0.6  0.24  0.97  0.23 0.56  0.46  0.61  0.61  0.36 2.1 9.5 10  19  6.5  18 1.2  6.2  26 0.006 1.2 3.5 0.62 7.2

(id="DEDGD,EMEES:DE.DEF.DEFG.DEFGK")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.T*2*c*T/(Kc*KcT)
            )*It
          +((kt.TT*2*C/KC)+2*kc.CC*C^2/(KC*KCC)+kt.TT*2*T*C/(KT*KTC) 
            + kt.TT*2*t*C/(Kt*KtC)
            + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
          + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
          )/2, 
    control=list(minFactor=1e-16),
    start=cd24[-c(3,5,6,8,9)],
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#kt.T kt.TT kt.TC kt.tT kt.cT kc.C kc.CC kc.TC kc.tC kc.cC KT   KTT KC  KCC KTC KtT KcT KtC  KcC  Kt  Kc   Ktt Ktc Kcc
#0.28   0.6  0.24  0.97  0.23 0.56  0.46  0.61  0.61  0.36 2.1 9.5 10  19  6.5  18 1.2  6.2  26 0.006 1.2 3.5 0.62 7.2

(id="DEDGD,LLEES:DE.DEF.DEFG.DEFGK")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KTC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.T*2*c*T/(Kc*KcT)
            )*It
          +((kc.C*2*C/KC)+2*kc.C*C^2/(KC*KCC)+kt.TT*2*T*C/(KT*KTC) 
            + kt.TT*2*t*C/(Kt*KtC)
            + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
          + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
          )/2, 
    control=list(minFactor=1e-16),
    start=cd24[-c(3,5,7,8,9)],
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#kt.T kt.TT kt.TC kt.tT kt.cT kc.C kc.CC kc.TC kc.tC kc.cC KT   KTT KC  KCC KTC KtT KcT KtC  KcC  Kt  Kc   Ktt Ktc Kcc
#0.28   0.6  0.24  0.97  0.23 0.56  0.46  0.61  0.61  0.36 2.1 9.5 10  19  6.5  18 1.2  6.2  26 0.006 1.2 3.5 0.62 7.2




(id="DDFDXFPPPX:DE.DID.DE.DEFI.III")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC) 
            +kt.tT*2*t*T/(Kt*KtT)
            +kt.tT*2*c*T/(Kc*KcT)
            )*It
          +((kt.T*2*C/KC)+kt.TC*2*T*C/(KT*KC) 
            + kt.tT*2*t*C/(Kt*KtC)
#            + kt.cT*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC + #C^2/(KC*KC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) # + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kc+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(
        kt.T=.88,
        kt.TC=.42,
        KT=12.1,
        KTT=144,
        KC=18,
#    
        kt.tT=0.4,
#        kt.cT=0.4,
        Kt=1.3,  
        Kc=2.7,
#        
        KtT=31,
        KcT=24,
        KtC=9
#        KcC=Inf
#        Kcc=Inf
#        Ktc=Inf
#        Ktt=Inf
    ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=TRUE)



(id="DDFDXFFFFX:DE.DID.DE.DEFI.III")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC) 
            +kt.TC*2*t*T/(Kt*KtT)
            +kt.TC*2*c*T/(Kc*KcT)
            )*It
          +((kt.T*2*C/KC)+kt.TC*2*T*C/(KT*KC) 
            + kt.TC*2*t*C/(Kt*KtC)
#            + kt.cT*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC + #C^2/(KC*KC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) # + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kc+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(
        kt.T=.88,
        kt.TC=.42,
        KT=12.1,
        KTT=144,
        KC=18,
#    
#        kt.tT=0.4,
#        kt.cT=0.4,
        Kt=1.3,  
        Kc=2.7,
#        
        KtT=31,
        KcT=24,
        KtC=9
#        KcC=Inf
#        Kcc=Inf
#        Ktc=Inf
#        Ktt=Inf
    ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=TRUE)


(id="DDFDXFFFFX:DE.DID.DD.DEFI.III")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC) 
            +kt.TC*2*t*T/(Kt*KtT)
            +kt.TC*2*c*T/(Kt*KcT)
            )*It
          +((kt.T*2*C/KC)+kt.TC*2*T*C/(KT*KC) 
            + kt.TC*2*t*C/(Kt*KtC)
#            + kt.cT*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC + #C^2/(KC*KC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kt*KcT) + 2*t*C/(Kt*KtC) # + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kt+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(
        kt.T=.88,
        kt.TC=.42,
        KT=12.1,
        KTT=144,
        KC=18,
#    
#        kt.tT=0.4,
#        kt.cT=0.4,
        Kt=2,  
#        Kc=2.7,
#        
        KtT=31,
        KcT=24,
        KtC=9
#        KcC=Inf
#        Kcc=Inf
#        Ktc=Inf
#        Ktt=Inf
    ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)


(id="DDFDXFFFFX:DE.DID.DD.DDFI.III")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC) 
            +kt.TC*2*t*T/(Kt*KtT)
            +kt.TC*2*c*T/(Kt*KtT)
            )*It
          +((kt.T*2*C/KC)+kt.TC*2*T*C/(KT*KC) 
            + kt.TC*2*t*C/(Kt*KtC)
#            + kt.cT*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC + #C^2/(KC*KC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kt*KtT) + 2*t*C/(Kt*KtC) # + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kt+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(
        kt.T=.88,
        kt.TC=.42,
        KT=12.1,
        KTT=144,
        KC=18,
#    
#        kt.tT=0.4,
#        kt.cT=0.4,
        Kt=2,  
#        Kc=2.7,
#        
        KtT=20,
#        KcT=20,
        KtC=9
#        KcC=Inf
#        Kcc=Inf
#        Ktc=Inf
#        Ktt=Inf
    ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DDFDXFFFFX:DE.SIS.DD.SSFI.III")
mod<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KC) 
            +kt.TC*2*t*T/(Kt*KC)
            +kt.TC*2*c*T/(Kt*KC)
            )*It
          +((kt.T*2*C/KC)+kt.TC*2*T*C/(KT*KC) 
            + kt.TC*2*t*C/(Kt*KtC)
#            + kt.cT*2*c*C/(Kc*KcC)
            )*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC + #C^2/(KC*KC) +
          2*t*T/(Kt*KC) + 2*c*T/(Kt*KC) + 2*t*C/(Kt*KtC) # + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kt+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(
        kt.T=.88,
        kt.TC=.42,
        KT=12.1,
        KTT=144,
        KC=18,
#    
#        kt.tT=0.4,
#        kt.cT=0.4,
        Kt=2,  
#        Kc=2.7,
#        
#        KtT=20,
#        KcT=20,
        KtC=9
#        KcC=Inf
#        Kcc=Inf
#        Ktc=Inf
#        Ktt=Inf
    ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)
names(tab)[1]<-"Model"

#Note:  I tried freeing up Ktc and Ktt at this point and these still made things singular

library(hwriter)
hwrite(t(tab),"C:/Users/radivot/case/active/papers/TK2/allDimer.html",cellspacing=0,cellpadding=2)


#
#plotmod<-function(mod,d,tab,id="none",jj=1,kill=FALSE,pause=FALSE){
#  if (class(mod)=="nls"){
#    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
#    if (DNA=="wt") {Cg=c(0,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
#    pd=rbind(data.frame(T=rep(Xg,length(Cg)),C=rep(Cg,each=length(Xg)),t=0,c=0,It=1,Ic=0),
#        data.frame(T=rep(Tg,each=length(Xg)),C=rep(Xg,length(Tg)),t=0,c=0,It=0,Ic=1) )
#    nX=length(Xg)
#    if (DNA=="wt") {
#      tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
##          tinh=list(list(t=0,c=c(1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0))
#      cinh=list(list(t=0,c=c(1,5)),list(t=5,c=0)) } else { 
##          tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
##      cinh=list(list(t=0,c=c(1,5)),list(t=1,c=0),list(t=5,c=0)) } else { 
#      tinh=list(list(t=0,c=c(1,10)),list(t=1,c=0),list(t=10,c=0))
#      cinh=list(list(t=0,c=c(2,10)),list(t=2,c=0),list(t=10,c=0)) 
#    }
#    for (i in 1:length(tinh)) {
#      nc=length(tinh[[i]]$c)
#      pd=rbind(pd,data.frame(T=rep(Xg,nc),C=0,t=tinh[[i]]$t,c=rep(tinh[[i]]$c,each=nX),It=1,Ic=0))
#    }
#    for (i in 1:length(cinh)) {
#      nc=length(cinh[[i]]$c)
#      pd=rbind(pd,data.frame(T=0,C=rep(Xg,nc),t=cinh[[i]]$t,c=rep(cinh[[i]]$c,each=nX),It=0,Ic=1))
#    }
##    pd=pd[-c(1:16,49:64),]
#    pd=data.frame(pd,Ek=predict(mod,pd))
#    if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
#    if (.Platform$OS.type=="windows") 
#      windows(width = 8, height = 8,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
#    par(mfrow=c(2,2),mar=c(5.2,4,0,1)+.1,oma=c(0,0,3,1))
#    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
#    for (i in 1:length(Cg)) 
#    {
#      with(d[(d$It>0)&(d$C==Cg[i])&(d$c==0)&(d$t==0),],points(T,k,col=i,pch=i))
#      with(pd[(pd$It>0)&(pd$C==Cg[i])&(pd$c==0)&(pd$t==0),],lines(T,Ek,col=i))
#    }
#    legend("topleft",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
#    
#    print(s<-summary(mod))
#    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
#    for (i in 1:length(Tg)) 
#    {
#      with(d[(d$Ic>0)&(d$T==Tg[i])&(d$c==0)&(d$t==0),],points(C,k,col=i,pch=i))
#      with(pd[(pd$Ic>0)&(pd$T==Tg[i])&(pd$c==0)&(pd$t==0),],lines(C,Ek,col=i))
#    }
#    legend("topleft",legend=paste("[dT]=",Tg,sep=""),pch=1:length(Tg),col=1:length(Tg),bty="n")
#    (p=dim(s$parameters)[1])
#    ed=data.frame(d,Ek=predict(mod,d))
#    SSE=sum((ed$Ek-ed$k)^2)
#    N=dim(d)[1]
#    P=p+1
#    sig=format(sqrt(SSE/(N-p)),digits=2)
#    AICc=format(N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N,digits=4)
#    SSE=format(SSE,digits=3)
##    AICuc=format(N*log(SSE/N)+2*P + N*log(2*pi) + N,digits=4)
##    AICchk=format(AIC(mod),digits=4)
#    title(paste(id," (",DNA,", AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=1.7,cex.main=1)
#    title(paste(sub(" ","",names(s$parameters[,1])),"=",format(s$parameters[,1],digits=1,trim=T,
#                nsmall=0,scientific=FALSE),sep="",collapse=";  "),
#        outer=TRUE,line=0.7,cex.main=0.8)
#    vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
#        "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
#    names(vals)<-names(s$parameters[,1])
#    tmp=c(id=sub(":DE.",":DE.\n",id),P=P,SSE=SSE,AIC=AICc,vals)
##tmp=c(id=id,p=p,aic=format(AIC(mod),digits=4),format(s$parameters[,1],digits=2))
#    if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
#    tab[jj,names(tmp)]=tmp
#    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
#    legs=NULL
#    for (i in 1:length(tinh)) 
#    {
#      for (j in 1:length(tinh[[i]]$c)) {
#        with(d[(d$It>0)&(d$t==tinh[[i]]$t)&(d$c==tinh[[i]]$c[j]),],points(T,k,col=i,pch=j))
#        with(pd[(pd$It>0)&(pd$t==tinh[[i]]$t)&(pd$c==tinh[[i]]$c[j]),],lines(T,Ek,col=i))
#        legs=rbind(legs,list(s=paste("[dTTP] =",format(tinh[[i]]$t,width=3),";  [dCTP] =",tinh[[i]]$c[j],sep=""),c=i,p=j))
#      }
#    }
#    legs=as.data.frame(legs)
#    legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
#    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase activity (1/sec)",xlab="[dC] uM",log="x")
#    legs=NULL
#    for (i in 1:length(cinh)) 
#    {
#      for (j in 1:length(cinh[[i]]$c)) {
#        with(d[(d$Ic>0)&(d$t==cinh[[i]]$t)&(d$c==cinh[[i]]$c[j]),],points(C,k,col=i,pch=j))
#        with(pd[(pd$Ic>0)&(pd$t==cinh[[i]]$t)&(pd$c==cinh[[i]]$c[j]),],lines(C,Ek,col=i))
#        legs=rbind(legs,list(s=paste("[dTTP] =",format(cinh[[i]]$t,width=3),";  [dCTP] =",cinh[[i]]$c[j],sep=""),c=i,p=j))
#      }
#    }
#    legs=as.data.frame(legs)
#    legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
#  } else {tab[jj,1]=id; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
#  if (pause) locator(n=1)
#  par(mfrow=c(1,1))
#  tab
#}

