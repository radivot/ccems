plotmod<-function(mod,tab,jj=1,kill=FALSE,pause=FALSE,showBad=FALSE,DNA="wt",
    Xg = c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)){
  id=mod$id
  d=mod$d
  if (class(mod)=="nls"){
    if (DNA=="wt") {if (!showBad) {Cg=c(0,10,20)} else Cg=c(0,5,10,20); Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
    pd=rbind(data.frame(T=rep(Xg,length(Cg)),C=rep(Cg,each=length(Xg)),It=1,Ic=0),
        data.frame(T=rep(Tg,each=length(Xg)),C=rep(Xg,length(Tg)),It=0,Ic=1) )
    pd=data.frame(pd,Ek=predict(mod,pd))
    if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
    if (.Platform$OS.type=="windows") 
      windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
    par(mfcol=c(1,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,3,1))
    plot(Xg,type="n",ylim=range(c(d$k)),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
    for (i in 1:length(Cg)) 
    {
      with(d[(d$It>0)&(d$C==Cg[i]),],points(T,k,col=i,pch=i))
      with(pd[(pd$It>0)&(pd$C==Cg[i]),],lines(T,Ek,col=i))
    }
    legend("topleft",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
    
    print(s<-summary(mod))
    plot(Xg,type="n",ylim=range(c(d$k)),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
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
    title(paste(id," (AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=1.7,cex.main=1)
    title(paste(sub(" ","",names(s$parameters[,1])),"=",format(s$parameters[,1],digits=1,trim=T),sep="",collapse=";  "),
        outer=TRUE,line=0.7,cex.main=0.8)
    if (pause) locator(n=1)
    par(mfrow=c(1,1))
    vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
        "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
#    vals=paste(format(s$parameters[,1],digits=2,trim=T),
#        "(",format(s$parameters[,2],digits=1,trim=T),")",sep="")
    names(vals)<-names(s$parameters[,1])
    tmp=c(id=id,p=p,SSE=SSE,AIC=AICc,vals)
#tmp=c(id=id,p=p,aic=format(AIC(mod),digits=4),format(s$parameters[,1],digits=2))
    if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
    tab[jj,names(tmp)]=tmp
  } else {tab[jj,1]=id; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  tab
}

plotmods<-function(mods,ds,picks,kill=FALSE,  Xg = c(.075,.1,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)){
  pd0=rbind(data.frame(T=Xg,C=0,It=1,Ic=0),
      data.frame(T=0,C=Xg,It=0,Ic=1) )
  if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
  if (.Platform$OS.type=="windows") 
    windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
  par(mfcol=c(1,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,0.5,0))
  plot(Xg,type="n",ylim=range(c(0,ds$k)),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
  for (j in 1:length(picks)) 
  {
    pd=data.frame(pd0,Ek=predict(mods[[picks[j]]],pd0))
    d=mods[[picks[j]]]$d
    with(d[d$It>0,],points(T,k,col=j,pch=j))
    with(pd[pd$It>0,],lines(T,Ek,col=j))
  }
  ids=sapply(mods,function(x) x$id)
  legend("topleft",legend=ids[picks], lty=1,col=1:length(picks),bty="n")
  plot(Xg,type="n",ylim=range(c(0,ds$k)),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
  for (j in 1:length(picks)) 
    if (class(mods[[picks[j]]])=="nls"){
      pd=data.frame(pd0,Ek=predict(mods[[picks[j]]],pd0))
      d=mods[[picks[j]]]$d
      with(d[d$Ic>0,],points(C,k,col=j,pch=j))
      with(pd[pd$Ic>0,],lines(C,Ek,col=j))
    }
  datasets=sapply(mods,function(x) x$data)[picks]
  legend("topleft",legend=datasets, pch=1:length(datasets),col=1:length(datasets),bty="n")
}

mkTab<-function(mods,ds,jj=1){
  options(stringsAsFactors = FALSE)
  (tab=data.frame(data=rep(0,length(mods)),id="XX",P=0,SSE=0,AIC=0,kt.T=0,kt.TT="kt.T",kc.C="kt.T",kc.CC="kt.T",KT=0,KTT=0,KC=0,KCC="Inf"))
  Xg <- c(.075,.1,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
  ids=sapply(mods,function(x) x$id)
  pd0=rbind(data.frame(T=Xg,C=0,It=1,Ic=0),
      data.frame(T=0,C=Xg,It=0,Ic=1) )
  for (j in 1:length(mods)) 
  {
    if (class(mods[[j]])=="nls"){
      pd=data.frame(pd0,Ek=predict(mods[[j]],pd0))
      d=mods[[j]]$d
#      d=subset(ds,year==mods[[j]]$year)
      s<-summary(mods[[j]])
      (p=dim(s$parameters)[1])
      ed=data.frame(d,Ek=predict(mods[[j]],d))
      SSE=sum((ed$Ek-ed$k)^2)
      N=dim(d)[1]
      P=p+1
      sig=format(sqrt(SSE/(N-p)),digits=2)
      AICc=format(N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N,digits=4)
      SSE=format(SSE,digits=3)
      vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
          "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
      names(vals)<-names(s$parameters[,1])
      tmp=c(data=mods[[j]]$data,id=ids[j],P=P,SSE=SSE,AIC=AICc,vals)
      if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
      tab[jj,names(tmp)]=tmp
    } else {tab[jj,1]=ids[j]; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  }
  tab
}

plotmodAll<-function(mod,d,tab,id="none",jj=1,od=NULL,kill=FALSE,pause=FALSE){
  if (class(mod)=="nls"){
    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
    if (DNA=="wt") {Cg=c(0,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
    pd=rbind(data.frame(T=rep(Xg,length(Cg)),C=rep(Cg,each=length(Xg)),t=0,c=0,It=1,Ic=0),
        data.frame(T=rep(Tg,each=length(Xg)),C=rep(Xg,length(Tg)),t=0,c=0,It=0,Ic=1) )
    nX=length(Xg)
    if (DNA=="wt") {
      tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
#          tinh=list(list(t=0,c=c(1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0))
      cinh=list(list(t=0,c=c(1,5)),list(t=5,c=0)) } else { 
#          tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
#      cinh=list(list(t=0,c=c(1,5)),list(t=1,c=0),list(t=5,c=0)) } else { 
      tinh=list(list(t=0,c=c(1,10)),list(t=1,c=0),list(t=10,c=0))
      cinh=list(list(t=0,c=c(2,10)),list(t=2,c=0),list(t=10,c=0)) 
    }
    for (i in 1:length(tinh)) {
      nc=length(tinh[[i]]$c)
      pd=rbind(pd,data.frame(T=rep(Xg,nc),C=0,t=tinh[[i]]$t,c=rep(tinh[[i]]$c,each=nX),It=1,Ic=0))
    }
    for (i in 1:length(cinh)) {
      nc=length(cinh[[i]]$c)
      pd=rbind(pd,data.frame(T=0,C=rep(Xg,nc),t=cinh[[i]]$t,c=rep(cinh[[i]]$c,each=nX),It=0,Ic=1))
    }
#    pd=pd[-c(1:16,49:64),]
    pd=data.frame(pd,Ek=predict(mod,pd))
    if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
    if (.Platform$OS.type=="windows") 
      windows(width = 9.5, height = 8,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
    par(mfrow=c(2,2),mar=c(5.5,4,0,1)+.1,oma=c(0,0,3,1))
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)", #xlab="[dT] uM",log="x")
        xlab=expression(paste("[dT] ",mu,"M")),log="x")
    with(od[(od$C==5),],points(T,k,col="orange",pch="X"))
    legend(.41,.48,"Deleted outliers at [dC]=5",pch="X",col="orange",bty="n")
    
    for (i in 1:length(Cg)) 
    {
      with(d[(d$It>0)&(d$C==Cg[i])&(d$c==0)&(d$t==0),],points(T,k,col=i,pch=i))
      with(pd[(pd$It>0)&(pd$C==Cg[i])&(pd$c==0)&(pd$t==0),],lines(T,Ek,col=i))
    }
    legend("topleft",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
    
    print(s<-summary(mod))
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",
        xlab=expression(paste("[dC] ",mu,"M")),log="x")
    for (i in 1:length(Tg)) 
    {
      with(d[(d$Ic>0)&(d$T==Tg[i])&(d$c==0)&(d$t==0),],points(C,k,col=i,pch=i))
      with(pd[(pd$Ic>0)&(pd$T==Tg[i])&(pd$c==0)&(pd$t==0),],lines(C,Ek,col=i))
    }
    legend("topleft",legend=paste("[dT]=",Tg,sep=""),pch=1:length(Tg),col=1:length(Tg),bty="n")
    (p=dim(s$parameters)[1])
    ed=data.frame(d,Ek=predict(mod,d))
    SSE=sum((ed$Ek-ed$k)^2)
    N=dim(d)[1]
    P=p+1
    sig=format(sqrt(SSE/(N-p)),digits=2)
    AICc=format(N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N,digits=4)
#    AICc=format(N*log(SSE/N)+2*P + N*log(2*pi) + N,digits=4) # checks out ok with AIC(mod)
    SSE=format(SSE,digits=3)
#    AICuc=format(N*log(SSE/N)+2*P + N*log(2*pi) + N,digits=4)
#    AICchk=format(AIC(mod),digits=4)
    title(paste(id," (AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=1.7,cex.main=1)
    title(paste(sub(" ","",names(s$parameters[,1])),"=",format(s$parameters[,1],digits=1,trim=T,
                nsmall=0,scientific=FALSE),sep="",collapse=";  "),
        outer=TRUE,line=0.7,cex.main=0.8)
    vals=paste(format(s$parameters[,1],digits=2,trim=T,scientific=FALSE),
        "(",format(s$parameters[,2],digits=1,trim=T,scientific=FALSE),")",sep="")
    names(vals)<-names(s$parameters[,1])
    tmp=c(id=id,p=p,SSE=SSE,AIC=AICc,vals)
#tmp=c(id=id,p=p,aic=format(AIC(mod),digits=4),format(s$parameters[,1],digits=2))
    if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
    tab[jj,names(tmp)]=tmp
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",#xlab="[dT] uM",log="x")
        xlab=expression(paste("[dT] ",mu,"M")),log="x")
    with(od[(od$C!=5),],points(T,k,col="orange",pch="X"))
    legend(.40,.5,"Deleted outliers at [dTTP]=1 and 10",pch="X",col="orange",bty="n")
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
    plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase activity (1/sec)", #xlab="[dC] uM",log="x")
        xlab=expression(paste("[dC] ",mu,"M")),log="x")
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
  } else {tab[jj,1]=id; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  if (pause) locator(n=1)
  par(mfrow=c(1,1))
  tab
}

plotFracA<-function(d){
  Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
  Cg=c(0,10,20);
  Tg=c(0,1,5)
  tinh=list(list(t=0,c=1),list(t=2,c=0),list(t=10,c=0))
  cinh=list(list(t=0,c=c(1,5)),list(t=5,c=0)) 
  windows(width = 9.5, height = 8,restoreConsole = TRUE,ypos=0)
  par(mfrow=c(2,2),mar=c(5.5,4,0,1)+.1,oma=c(0,0,0,0))
  plot(Xg,type="n",ylim=range(d$frcA),xlim=range(Xg),ylab="dT Kinase activity (proportion)", #xlab="[dT] uM",log="x")
      xlab=expression(paste("[dT] ",mu,"M")),log="x")
  
  for (i in 1:length(Cg)) 
    with(d[(d$It>0)&(d$C==Cg[i])&(d$c==0)&(d$t==0),],points(T,frcA,col=i,pch=i))
  legend("topright",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
  plot(Xg,type="n",ylim=range(d$frcA),xlim=range(Xg),ylab="dC Kinase activity (proportion)",
      xlab=expression(paste("[dC] ",mu,"M")),log="x")
  for (i in 1:length(Tg)) 
    with(d[(d$Ic>0)&(d$T==Tg[i])&(d$c==0)&(d$t==0),],points(C,frcA,col=i,pch=i))
  legend("topright",legend=paste("[dT]=",Tg,sep=""),pch=1:length(Tg),col=1:length(Tg),bty="n")
  
  plot(Xg,type="n",ylim=range(d$frcA),xlim=range(Xg),ylab="dT Kinase activity (proportion)",#xlab="[dT] uM",log="x")
      xlab=expression(paste("[dT] ",mu,"M")),log="x")
  legs=NULL
  for (i in 1:length(tinh)) 
  {
    for (j in 1:length(tinh[[i]]$c)) {
      with(d[(d$It>0)&(d$t==tinh[[i]]$t)&(d$c==tinh[[i]]$c[j]),],points(T,frcA,col=i,pch=j))
      legs=rbind(legs,list(s=paste("[dTTP] =",format(tinh[[i]]$t,width=3),";  [dCTP] =",tinh[[i]]$c[j],sep=""),c=i,p=j))
    }
  }
  legs=as.data.frame(legs)
  legend("topright",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
  plot(Xg,type="n",ylim=range(d$frcA),xlim=range(Xg),ylab="dC Kinase activity (proportion)", #xlab="[dC] uM",log="x")
      xlab=expression(paste("[dC] ",mu,"M")),log="x")
  legs=NULL
  for (i in 1:length(cinh)) 
  {
    for (j in 1:length(cinh[[i]]$c)) {
      with(d[(d$Ic>0)&(d$t==cinh[[i]]$t)&(d$c==cinh[[i]]$c[j]),],points(C,frcA,col=i,pch=j))
      legs=rbind(legs,list(s=paste("[dTTP] =",format(cinh[[i]]$t,width=3),";  [dCTP] =",cinh[[i]]$c[j],sep=""),c=i,p=j))
    }
  }
  legs=as.data.frame(legs)
  legend("topright",legend=unlist(legs$s),pch=unlist(legs$p),col=unlist(legs$c),cex=0.9,bty="n")
  par(mfrow=c(1,1))
}




ARD12<-function(ardP)
{ attach(ardP)
  KT =2/(1/KAT + 1/KBT)
  KC =2/(1/KAC + 1/KBC)
  KCC=KAC*KBC/KC
  KTT=KAT*KBT/KT
  KTC=2/(1/(KAC*KBT) + 1/(KAT*KBC))/KT
  Kc =2/(1/KAc + 1/KBc)
  Kt =2/(1/KAt + 1/KBt)
  Kcc=KAc*KBc/Kc
  Ktt=KAt*KBt/Kt
  Ktc=2/(1/(KAc*KBt) + 1/(KAt*KBc))/Kt
  KtT=2/(1/(KAT*KBt) + 1/(KAt*KBT))/Kt
  KtC=2/(1/(KAC*KBt) + 1/(KAt*KBC))/Kt
  KcT=2/(1/(KAc*KBT) + 1/(KAT*KBc))/Kc
  KcC=2/(1/(KAc*KBC) + 1/(KAC*KBc))/Kc
  kt.T=(KT/2)*(kt.A/KAT + kt.B/KBT)
  kt.TT=(KTT*KT/2)*(kt.A+kt.B)/(KAT*KBT)
  kt.TC=(KTC*KT/2)*(kt.A/(KAT*KBC) + kt.B/(KAC*KBT))
  kt.tT=(KtT*Kt/2)*(kt.A/(KAT*KBt) + kt.B/(KAt*KBT))
  kt.cT=(KcT*Kc/2)*(kt.A/(KAT*KBc) + kt.B/(KAc*KBT))
  kc.C=(KC/2)*(kc.A/KAC + kc.B/KBC)
  kc.CC=(KCC*KC/2)*(kc.A+kc.B)/(KAC*KBC)
  kc.TC=(KTC*KT/2)*(kc.A/(KAC*KBT) + kc.B/(KBC*KAT))
  kc.tC=(KtC*Kt/2)*(kc.A/(KAC*KBt) + kc.B/(KAt*KBC))
  kc.cC=(KcC*Kc/2)*(kc.A/(KAC*KBc) + kc.B/(KAc*KBC))
  detach(ardP)
  list( kt.T =kt.T, kt.TT=kt.TT, kt.TC=kt.TC, kt.tT=kt.tT, kt.cT=kt.cT, 
      kc.C =kc.C, kc.CC=kc.CC, kc.TC=kc.TC, kc.tC=kc.tC, kc.cC=kc.cC,
      KT =KT, KTT=KTT,  KtT=KtT, KcT=KcT, KC =KC, KCC=KCC, KTC=KTC, KtC=KtC,
      KcC=KcC, Kt =Kt, Kc =Kc, Ktt=Ktt, Ktc=Ktc, Kcc=Kcc )
}

ARD8<-function(ardP)
{ attach(ardP)
  KT =2/(1/KAT + 1/KBT)
  KC =2/(1/KAC + 1/KBC)
  KCC=KAC*KBC/KC
  KTT=KAT*KBT/KT
  KTC=2/(1/(KAC*KBT) + 1/(KAT*KBC))/KT
  kt.T=(KT/2)*(kt.A/KAT + kt.B/KBT)
  kt.TT=(KTT*KT/2)*(kt.A+kt.B)/(KAT*KBT)
  kt.TC=(KTC*KT/2)*(kt.A/(KAT*KBC) + kt.B/(KAC*KBT))
  kc.C=(KC/2)*(kc.A/KAC + kc.B/KBC)
  kc.CC=(KCC*KC/2)*(kc.A+kc.B)/(KAC*KBC)
  kc.TC=(KTC*KT/2)*(kc.A/(KAC*KBT) + kc.B/(KBC*KAT))
  detach(ardP)
  list( kt.T =kt.T, kt.TT=kt.TT, kt.TC=kt.TC, 
      kc.C =kc.C, kc.CC=kc.CC, kc.TC=kc.TC,
      KT =KT, KTT=KTT, KC =KC, KCC=KCC, KTC=KTC)
}

getEk<-function(g,init=FALSE){ 
  if (sum(g$params[,"constr"]!="none")>0) {
    bound=row.names(g$params)[g$params[,"constr"]!="none"]  # bound = follower, free = leader
    free=g$params[bound,"constr"] # these are the free params to which the bound ones are bound
    g$params[bound,"final"]=g$params[free,"final"]
  }
  Pframe=as.data.frame(t(g$params[,"final",drop=F]))
  if (g$CD) {
    with(Pframe,g$d<<-transform(g$d,  # borrow CD for Eq.2
            Ek= ((kt.A*(T/KAT)^h/(1+(T/KAT)^h+t/KAt+c/KAc)/(1+C/KAC)))*It
                +((kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)/(1+t/KBt)))*Ic
        )  ) } else {
#            Ek=
#                ((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT)+kt.TC*2*T*C/(KT*KTC) 
#                    +kt.tT*2*t*T/(Kt*KtT)
#                    +kt.cT*2*c*T/(Kc*KcT)
#                    )*It
#                  +((kc.C*2*C/KC)+2*kc.CC*C^2/(KC*KCC)+kc.TC*2*T*C/(KT*KTC) 
#                    + kc.tC*2*t*C/(Kt*KtC)
#                    + kc.cC*2*c*C/(Kc*KcC)
#                    )*Ic)/
#                (1+2*T/KT+T^2/(KT*KTT)+2*T*C/(KT*KTC)+2*C/KC + C^2/(KC*KCC) +
#                  2*t*T/(Kt*KtT) + 2*c*T/(Kc*KcT) + 2*t*C/(Kt*KtC) + 2*c*C/(Kc*KcC) # terms from numerator that may have activity
#                  + 2*c/Kc+2*t/Kt  + t*t/(Kt*Ktt) + 2*t*c/(Kt*Ktc) + c*c/(Kc*Kcc)
#                  )/2
#        )  ) } else {
#    with(Pframe,g$d<<-transform(g$d,
#            Ek=((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#                  +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It/2
#                +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#                  +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2
#        )  )
    with(Pframe,g$d<<-transform(g$d,
            Ek=((kt.Af*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
                  +(kt.Bf*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It
                +((kc.Af*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
                  +(kc.Bf*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic
        )  )
    with(Pframe,g$d<<-transform(g$d,
            EkA=((kt.Af*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc)) )*It
                +((kc.Af*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc)))*Ic
        )  )
    with(Pframe,g$d<<-transform(g$d,
            EkB=((kt.Bf*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It
                +((kc.Bf*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic
        )  )
  }
  SSE=sum((g$d$Ek-g$d$k)^2)
  N=g$nData=dim(g$d)[1]
  P=sum(g$params[,"opt"]) + 1 # include the variance
  aic=N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N
  if (P>=N-1) aic=Inf
  if (init) {
    g$SSE$initial=SSE 
    g$AIC$initial = aic
  } 
  g$SSE$final=SSE 
  g$AIC$final =  aic
  g
}


fitModel <- function(model) {
  
  fopt <- function(pars,model) {
    pars=exp(pars)
    model$params[names(pars),"final"]=pars  
    model=getEk(model)
    return(model$SSE$final)   
  }
  
  p0=t(model$params[model$params[,"opt"],"initial",drop=F])
  p0=unlist(as.data.frame(p0))
  if ((model$nOptParams<-length(p0))>0) {                    
    p0=sapply(p0,log)
    if (model$nOptParams>1) opt<-optim(p0,fopt,hessian=TRUE,model=model,control=list(trace=F,maxit=5000)) else
      opt<-optim(p0,fopt,method="BFGS",hessian=TRUE,model=model,control=list(trace=F));       
  }
  opar<-exp(opt$par)
  model$params[names(opar),"final"]=opar  
  model=getEk(model)
  sg<-sqrt(opt$value/(model$nData-model$nOptParams))
  if (det(opt$hessian)>0) 
  {sig=sg*sqrt(diag(solve(opt$hessian/2)))
    upper=signif(opt$par+1.96*sig,3)
    lower=signif(opt$par-1.96*sig,3)
    CI=cbind(lower,upper)
    model$CI=exp(CI); 
    model$hess=TRUE} else
  {model$hess=FALSE} 
  model
#  opt
}




plotAllg<-function(g,tab,jj=1,kill=FALSE){
  Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
  DNA=g$DNA
  id=g$id
  if (DNA=="wt") {Cg=c(0,10,20);Tg=c(0,1,5)} else {Cg=c(0,1,10);Tg=c(0,2,10)}
  pd=rbind(data.frame(T=rep(Xg,length(Cg)),C=rep(Cg,each=length(Xg)),t=0,c=0,It=1,Ic=0),
      data.frame(T=rep(Tg,each=length(Xg)),C=rep(Xg,length(Tg)),t=0,c=0,It=0,Ic=1) )
  nX=length(Xg)
  if (DNA=="wt") {
    tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
#          tinh=list(list(t=0,c=c(1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0))
    cinh=list(list(t=0,c=c(1,5)),list(t=5,c=0)) } else { 
#          tinh=list(list(t=0,c=c(1)),list(t=2,c=0),list(t=10,c=0))
#      cinh=list(list(t=0,c=c(1,5)),list(t=1,c=0),list(t=5,c=0)) } else { 
    tinh=list(list(t=0,c=c(1,10)),list(t=1,c=0),list(t=10,c=0))
    cinh=list(list(t=0,c=c(2,10)),list(t=2,c=0),list(t=10,c=0)) 
  }
  for (i in 1:length(tinh)) {
    nc=length(tinh[[i]]$c)
    pd=rbind(pd,data.frame(T=rep(Xg,nc),C=0,t=tinh[[i]]$t,c=rep(tinh[[i]]$c,each=nX),It=1,Ic=0))
  }
  for (i in 1:length(cinh)) {
    nc=length(cinh[[i]]$c)
    pd=rbind(pd,data.frame(T=0,C=rep(Xg,nc),t=cinh[[i]]$t,c=rep(cinh[[i]]$c,each=nX),It=0,Ic=1))
  }
  pg=g
  pg$d=pd
  pg=getEk(pg)
  pd=pg$d
  d=g$d
  if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
  if (.Platform$OS.type=="windows") 
    windows(width = 9.5, height = 8,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
  par(mfrow=c(2,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,4,1))
  plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
  for (i in 1:length(Cg)) 
  {
    with(d[(d$It>0)&(d$C==Cg[i])&(d$c==0)&(d$t==0),],points(T,k,col=i,pch=i))
    with(pd[(pd$It>0)&(pd$C==Cg[i])&(pd$c==0)&(pd$t==0),],lines(T,Ek,col=i))
  }
  legend("topleft",legend=paste("[dC]=",Cg,sep=""),pch=1:length(Cg),col=1:length(Cg),bty="n")
  plot(Xg,type="n",ylim=range(d$k),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
  for (i in 1:length(Tg)) 
  {
    with(d[(d$Ic>0)&(d$T==Tg[i])&(d$c==0)&(d$t==0),],points(C,k,col=i,pch=i))
    with(pd[(pd$Ic>0)&(pd$T==Tg[i])&(pd$c==0)&(pd$t==0),],lines(C,Ek,col=i))
  }
  legend("topleft",legend=paste("[dT]=",Tg,sep=""),pch=1:length(Tg),col=1:length(Tg),bty="n")
  SSE=g$SSE$final
  AICc=signif(g$AIC$final,3) 
  P=sum(g$params[,"opt"]) + 1 # include the variance
  SSE=format(SSE,digits=3)
  title(paste(id," (",DNA,", AIC=",AICc,", SSE=",SSE,")",sep=""),outer=TRUE,line=2.7,cex.main=1)
  title(paste(row.names(g$params)[1:10],"=",format(g$params[1:10,"final"],digits=2,trim=T,
              nsmall=0,scientific=FALSE),sep="",collapse=";  "), outer=TRUE,line=1.7,cex.main=0.8)
  title(paste(row.names(g$params)[11:24],"=",format(g$params[11:24,"final"],digits=2,trim=T,
              nsmall=0,scientific=FALSE),sep="",collapse=";  "), outer=TRUE,line=0.7,cex.main=0.8)
  vals=format(g$params[,"final"],digits=2,trim=T,scientific=FALSE)
#    vals=paste(format(g$params[,"final"],digits=2,trim=T,scientific=FALSE),
#        format(g$params[,"CI"],digits=1,trim=T,scientific=FALSE),sep="")
  names(vals)<-row.names(g$params)
  tmp=c(id=id,p=P-1,SSE=SSE,AIC=AICc,vals)
  if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
  tab[jj,names(tmp)]=tmp
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
  par(mfrow=c(1,1))
  tab
}
