rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,ET=E/2)
#TK2
DNA="H121N"
DNA="wt"
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(frstAut=="Wang")&(seq==DNA),select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[2:3]=c("T","C")
d=transform(d,It=as.numeric(!is.na(d$VdT)),Ic=as.numeric(!is.na(d$VdC)))
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
#d

plotmod<-function(mod,d,tab,id,jj=1,kill=FALSE,pause=FALSE){
  if (class(mod)=="nls"){
    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
    if (DNA=="wt") {Cg=c(0,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
#    if (DNA=="wt") {Cg=c(0,5,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
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
    vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
        "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
#    vals=paste(format(s$parameters[,1],digits=2,trim=T),
#        "(",format(s$parameters[,2],digits=1,trim=T),")",sep="")
    names(vals)<-names(s$parameters[,1])
    tmp=c(id=id,P=P,SSE=SSE,AIC=AICc,vals)
#tmp=c(id=id,p=p,aic=format(AIC(mod),digits=4),format(s$parameters[,1],digits=2))
    if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
    tab[jj,names(tmp)]=tmp
    } else {tab[jj,1]=id; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  tab
}

options(stringsAsFactors = FALSE)
############# START WITH 6 IN TABLE BY LOGICAL PROGRESSION ###############

# ONE MONOMER
(tab=data.frame(id="XX",p=rep(0,1),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kc.C="kt",KT=0,KC="KT"))
(id="DE:DE")
mod<-nls(k~(kt.T*T/KT/(1+T/KT+C/KC))*It
          +(kc.C*C/KC/(1+T/KT+C/KC))*Ic, 
    start=list(
        kt.T=.6,
        kc.C=.6,
        KT=5,
        KC=7),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=T)
tab1=tab
tab1
library(hwriter)
hwrite(tab1,"C:/Users/radivot/case/active/papers/TK2/oneMonomer.html")

# Asymmetric Rigid Dimers (ARDs)
(tab=data.frame(id="XX",P=0,SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.A=0,kt.B="kt.A",kc.A="X",kc.B=c(rep("kc.A",3)),
          KAT=0,KBT=0,KAC="Inf",KBC="KAC"))

(id="DEFG:DE.DE")
mod<-nls(k~(kt.A*T/KAT/(1+T/KAT+C/KAC) + kt.B*T/KBT/(1+T/KBT+C/KBC) )*It/2
          +(kc.A*C/KAC/(1+T/KAT+C/KAC) + kc.B*C/KBC/(1+T/KBT+C/KBC) )*Ic/2, 
    start=list(
        kt.A=.5,kt.B=.8,
        kc.A=.5,kc.B=.3,
        KAT=3,
        KBT=46,
        KAC=13,
        KBC=5),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEFF:DE.DD")
mod<-nls(k~(kt.A*T/KAT/(1+T/KAT+C/KAC) + kt.B*T/KBT/(1+T/KBT+C/KAC) )*It/2
        +(kc.A*C/KAC/(1+T/KAT+C/KAC) + kc.A*C/KAC/(1+T/KBT+C/KAC) )*Ic/2, 
    start=list(
        kt.A=.6,kt.B=1,
        kc.A=.5,
        KAT=3,KBT=100,
        KAC=8),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)

(id="DEFX:DE.ID")
mod<-nls(k~(kt.A*T/KAT/(1+T/KAT) + kt.B*T/KBT/(1+T/KBT+C/KBC) )*It/2
        +(                         kc.B*C/KBC/(1+T/KBT+C/KBC) )*Ic/2, 
    start=list(
        kt.A=.6,kt.B=1,
        kc.B=.5,
        KAT=3,KBT=100,
        KBC=8),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id)


tab1=tab
tab1
names(tab1)[1]<-"Model"
#names(tab1)[9:12]<-c("KA_T","KB_T","KA_C","KB_C")
library(hwriter)
hwrite(tab1,"C:/Users/radivot/case/active/papers/TK2/ARDTC.html",cellspacing=0,cellpadding=2)



#
#(id="DDFF:DE.DE:p")
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.A*T/KBT/(1+T/KBT+C/KBC)))*It
#        +(p*(kc.A*C/KAC/(1+T/KAT+C/KAC))+(1-p)*(kc.A*C/KBC/(1+T/KBT+C/KBC)))*Ic, 
#    start=list(
#        kt.A=.6,
#        kc.A=.6,
#        KAT=3,KBT=100,
#        KAC=8,KBC=8,p=0.5),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)

#(id="DEDD:DE.DE:p")
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.B*T/KBT/(1+T/KBT+C/KBC)))*It
#        +(p*(kt.A*C/KAC/(1+T/KAT+C/KAC))+(1-p)*(kt.A*C/KBC/(1+T/KBT+C/KBC)))*Ic, 
#    start=list(
#        kt.A=.5,kt.B=.8,
#        KAT=3,KBT=46,
#        KAC=13,KBC=5,p=0.4),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#
#(id="DEDD:DE.DD:p")
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.B*T/KBT/(1+T/KBT+C/KAC)))*It
#        +(p*(kt.A*C/KAC/(1+T/KAT+C/KAC))+(1-p)*(kt.A*C/KAC/(1+T/KBT+C/KAC)))*Ic, 
#    start=list(
#        kt.A=.5,kt.B=.8,
#        KAT=3,KBT=46,
#        KAC=13,p=0.5),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#
#(id="DDDD:DE.DE:p")
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.A*T/KBT/(1+T/KBT+C/KBC)))*It
#        +(p*(kt.A*C/KAC/(1+T/KAT+C/KAC))+(1-p)*(kt.A*C/KBC/(1+T/KBT+C/KBC)))*Ic, 
#    start=list(
#        kt.A=.6,
#        KAT=3,KBT=100,
#        KAC=8,KBC=8,p=0.5),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)
#
#(id="DDDX:DE.DI:p") # 
##p=.5
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.A*T/KBT/(1+T/KBT)))*It
#        +(p*(kt.A*C/KAC/(1+T/KAT+C/KAC)))*Ic, 
#    start=list(
#        kt.A=.5, #kt.B=.8,
##        kc.A=.5,
#        KAT=3,KBT=46,
#        KAC=8,p=.5),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id)


# THESE FAILED TO CONVERGE
#(id="DEFX:DEFI")
#p=.5
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.B*T/KBT/(1+T/KBT)))*It
#          +(p*(kc.A*C/KAC/(1+T/KAT+C/KAC)))*Ic, 
#    start=list(
#        kt.A=.6,kt.B=1,
#        kc.A=.5,
#        KAT=3,KBT=100,
#        KAC=8),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id,kill=T)
#
#(id="DEDX:DEFI")
#p=.5
#mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC))+(1-p)*(kt.B*T/KBT/(1+T/KBT)))*It
#          +(p*(kt.A*C/KAC/(1+T/KAT+C/KAC)))*Ic, 
#    start=list(
#        kt.A=.6,kt.B=1,
#        KAT=5,KBT=300,
#        KAC=7),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmod(mod,d,tab,id,kill=T)
