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
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(year==2003)&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
dtm=subset(d,(It==1)&(C==0))
dcm=subset(d,(Ic==1)&(T==0))
d=rbind(dtm,dcm)
d
killBad=FALSE
if (killBad) 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d



plotmod<-function(mod,d,tab,id="none",jj=1,kill=FALSE,pause=FALSE){
  if (class(mod)=="nls"){
    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
    nX=length(Xg)
    if (DNA=="wt") {
      if (!killBad) {
            tinh=list(list(t=0,c=c(0,1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0)) } else { 
            tinh=list(list(t=0,c=c(0,1)),list(t=2,c=0),list(t=10,c=0))}
      cinh=list(list(t=0,c=c(0,1,5)),list(t=1,c=0),list(t=5,c=0))  } else { 
      tinh=list(list(t=0,c=c(0,1,10)),list(t=1,c=0),list(t=10,c=0))
      cinh=list(list(t=0,c=c(0,2,10)),list(t=2,c=0),list(t=10,c=0)) 
    }

    pd=NULL
    for (i in 1:length(tinh)) {
      nc=length(tinh[[i]]$c)
      pd=rbind(pd,data.frame(T=rep(Xg,nc),C=0,t=tinh[[i]]$t,c=rep(tinh[[i]]$c,each=nX),It=1,Ic=0))
    }
    for (i in 1:length(cinh)) {
      nc=length(cinh[[i]]$c)
      pd=rbind(pd,data.frame(T=0,C=rep(Xg,nc),t=cinh[[i]]$t,c=rep(cinh[[i]]$c,each=nX),It=0,Ic=1))
    }
    pd=data.frame(pd,Ek=predict(mod,pd))
    if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
    if (.Platform$OS.type=="windows") 
      windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
    par(mfcol=c(1,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,3.5,1))
    
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
    legs=NULL
    for (i in 1:length(tinh)) 
    {
      for (j in 1:length(tinh[[i]]$c)) {
        with(d[(d$It>0)&(d$t==tinh[[i]]$t)&(d$c==tinh[[i]]$c[j]),],points(T,k,col=i,pch=j))
        with(pd[(pd$It>0)&(pd$t==tinh[[i]]$t)&(pd$c==tinh[[i]]$c[j]),],lines(T,Ek,col=i))
        legs=rbind(legs,list(s=paste("[dTTP] =",format(tinh[[i]]$t,width=3),";  [dCTP] =",tinh[[i]]$c[j],sep=""),c=i,p=j))
      }
    }
    legs=as.data.frame(legs)
    legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase activity (1/sec)",xlab="[dC] uM",log="x")
    legs=NULL
    for (i in 1:length(cinh)) 
    {
      for (j in 1:length(cinh[[i]]$c)) {
        with(d[(d$Ic>0)&(d$t==cinh[[i]]$t)&(d$c==cinh[[i]]$c[j]),],points(C,k,col=i,pch=j))
        with(pd[(pd$Ic>0)&(pd$t==cinh[[i]]$t)&(pd$c==cinh[[i]]$c[j]),],lines(C,Ek,col=i))
        legs=rbind(legs,list(s=paste("[dTTP] =",format(cinh[[i]]$t,width=3),";  [dCTP] =",cinh[[i]]$c[j],sep=""),c=i,p=j))
      }
    }
    legs=as.data.frame(legs)
    legend("topleft",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
    print(s<-summary(mod))
    (p=dim(s$parameters)[1])
    ed=data.frame(d,Ek=predict(mod,d))
    SSE=sum((ed$Ek-ed$k)^2)
#    print(ed)
    N=dim(d)[1]
    P=p+1
    sig=format(sqrt(SSE/(N-p)),digits=2)
    AICc=format(N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N,digits=4)
    SSE=format(SSE,digits=3)
#    AICuc=format(N*log(SSE/N)+2*P + N*log(2*pi) + N,digits=4)
#    AICchk=format(AIC(mod),digits=4)
    title(paste(id," (",DNA,", AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=2.2,cex.main=1)
    title(paste(sub(" ","",names(s$parameters[,1])),"=",format(s$parameters[,1],digits=1,trim=T),sep="",collapse="   "),
        outer=TRUE,line=1.2,cex.main=0.8)
    if (pause) locator(n=1)
    par(mfrow=c(1,1))
    vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
        "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
#    vals=paste(format(s$parameters[,1],digits=2,trim=T),"(",format(s$parameters[,2],digits=1,trim=T),")",sep="")
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
(tab=data.frame(id="XX",p=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.tT=0,kt.cT=0,
          kc.cC="kt.tT",kc.tC=0,
          KtT=0,KcT=0,
          KtC=0,KcC="KtT",
          Kt=0,Ktt=Inf,
          Kc="Kt",Kcc=Inf))


#0.88(0.02)  kt.T  0.42(0.04)  kt.T  X kt.TC 12.13(0.76) 144.34(18.17) 18.18(1.27) Inf KC
#DDFDXF:DE.DID 0.89 kt.T 0.40 kt.T X kt.TC 11.6 115 19 Inf KC 
kt.T=.88
kt.TT=kt.T
kt.TC=.42
kc.C=kt.T
kc.CC=0 # X actually
kc.TC=kt.TC
KT=12.1
KTT=144
KC=18
KCC=Inf
KTC=KC


(id="E0.E0:DE.DE.DI.DI")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC) +
            kt.tT*2*t*T/(Kt*KtT))*It
          +(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC)+ 
            kt.tT*2*t*C/(Kt*KtC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC) +
          2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          +2*c/Kc+2*t/Kt
          )/2, 
    control=list(minFactor=1e-16),
    start=list(kt.tT=0.4,
        KtT=19,KcT=4.3,
        KtC=3.1,KcC=7.6,
        Kt=1.3,
        Kc=2.7),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)



g=KC/KT

(id="HLLH:DE.DE.DEF.pE")
kt.tT=0.88
kt.cT=.42
kc.tC=.42
kt.cT=.0
kc.tC=.0
kc.cC=.88
Kcc=Inf
KcC=Inf
KcT=Inf
Ktt=Inf
Ktc=Inf
KtT=Inf
KtC=Inf
Kc=Inf
Kt=Inf




mod<-nls(k~((kt.T*2*T/KT + 2*kt.TT*T^2/(KT*KTT) + kt.TC*2*T*C/(KT*KC) 
#            + kt.tT*2*t*T/(Kt*KtT) + kt.cT*2*c*T/(Kc*KcT)
            )*It
          +(kc.C*2*C/KC + kc.TC*2*T*C/(KT*KC) #+2*kc.CC*C^2/(KC*KCC) 
#            + kc.tC*2*t*C/(Kt*KtC) + kc.cC*2*c*C/(Kc*KcC)
            )*Ic)/
        (1 + 2*T/KT + T^2/(KT*KTT) + 2*T*C/(KT*KC) #+2*C/KC+C^2/(KC*KCC) 
#          + 2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
          + 2*c/(Kc) + 2*t/Kt# +t^2/(Kt*Ktt)  + 2*t*c/(Kt*Ktc)+ c^2/(Kc*Kcc)
          )/2, 
    

#mod<-nls(k~((kt.T*2*T/KT + 2*kt.TT*T^2/(KT*KTT) + kt.TC*2*T*C/(KT*KC) +
#            kt.tT*2*t*T/(Kt*KtT) + kt.cT*2*c*T/(Kc*KcT))*It
#          +(kc.C*2*C/KC + kc.TC*2*T*C/(KT*KC) #+2*kc.CC*C^2/(KC*KCC) 
#            + kc.tC*2*t*C/(Kt*KtC) + kc.cC*2*c*C/(Kc*KcC))*Ic)/
#        (1 + 2*T/KT + T^2/(KT*KTT) + 2*T*C/(KT*KC) #+2*C/KC+C^2/(KC*KCC) 
#          + 2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
#          + 2*c/(Kc) + 2*t/Kt+ t^2/(Kt*Ktt)  + 2*t*c/(Kt*Ktc)+ c^2/(Kc*Kcc)
#          )/2, 
    
#mod<-nls(k~((kt.T*2*T/KT + 2*kt.TT*T^2/(KT*KTT) + kt.TC*2*T*C/(KT*KC) +
#            kt.tT*2*t*T/(Kt*KtT) + kt.cT*2*c*T/(g*Kt*KcT))*It
#          +(kc.C*2*C/KC + kc.TC*2*T*C/(KT*KC) #+2*kc.CC*C^2/(KC*KCC) 
#            + kc.tC*2*t*C/(Kt*KtC) + kc.cC*2*c*C/(g*Kt*KcC))*Ic)/
#        (1 + 2*T/KT + T^2/(KT*KTT) + 2*T*C/(KT*KC) #+2*C/KC+C^2/(KC*KCC) 
#          + 2*t*T/(Kt*KtT) + 2*c*T/(g*Kt*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(g*Kt*KcC) # terms from numerator: may have activity
#          + 2*c/(g*Kt) + 2*t/Kt+ t^2/(Kt*Ktt)  + 2*t*c/(Kt*Ktc)+ c^2/(g*Kt*Kcc)
#          )/2, 
    control=list(minFactor=1e-16),model=TRUE,
    start=list(
#        KtT=19,
#        KcT=3,
#        KtC=10,
#        Kc=100,
#        Kt=100
#        Ktt=10,
#        Ktc=1
        ),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=TRUE)






(id="D0.D0:DE.DE.DE.DI")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC) +
         kt.tT*2*t*T/(Kt*KtT))*It
        +(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC)+ 
        kc.tC*2*t*C/(Kt*KtC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC) +
        2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
        +2*c/Kc+2*t/Kt+t^2/(Kt*Ktt)
        )/2, 
    control=list(minFactor=1e-16),
    start=list(kt.tT=0.8,
               kc.tC=.4,
               KtT=11,KcT=3,
               KtC=10,KcC=1,
               Kt=1,Ktt=1,
               Kc=1),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="D0.D0:DE.DE.DI.DI")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC) +
         kt.tT*2*t*T/(Kt*KtT))*It
        +(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC)+ 
        kc.tC*2*t*C/(Kt*KtC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC) +
        2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator: may have activity
        +2*c/Kc+2*t/Kt
        )/2, 
    control=list(minFactor=1e-16),
    start=list(kt.tT=0.8,
               kc.tC=.4,
               KtT=11,KcT=3,
               KtC=10,KcC=1,
               Kt=1,
               Kc=1),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

#(id="D0.D0:DE.DE.DD.DI")  skip ... params want negative



(id="E0.E0:DE.DE.EI.EI")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC) +
         kt.tT*2*t*T/(Kt*KtT))*It
        +(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC)+ 
        kt.tT*2*t*C/(Kt*KtC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC) +
        2*t*T/(Kt*KtT) + 2*c*T/(Kt*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kt*KcC) # terms from numerator: may have activity
        +2*c/Kt+2*t/Kt
        )/2, 
    control=list(minFactor=1e-16),
    start=list(kt.tT=0.4,
               KtT=11,KcT=3,
               KtC=10,KcC=1,
               Kt=1),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="E0.E0:SE.DS.EI.EI")
mod<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.T*2*T*C/(KT*KC) +
         kt.tT*2*t*T/(Kt*KtT))*It
        +(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC)+kc.C*2*T*C/(KT*KC)+ 
        kt.tT*2*t*C/(Kt*KtC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KC)+2*C/KC+C^2/(KC*KC) +
        2*t*T/(Kt*KtT) + 2*c*T/(Kt*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kt*KtT) # terms from numerator: may have activity
        +2*c/Kt+2*t/Kt
        )/2, 
    control=list(minFactor=1e-16),
    start=list(kt.tT=0.4,
               KtT=11,KcT=3,
               KtC=10,
               Kt=1),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

library(hwriter)
hwrite(tab,"C:/Users/radivot/case/active/papers/TK2/tab3.html")

