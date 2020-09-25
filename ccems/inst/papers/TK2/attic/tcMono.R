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

plotmod<-function(mod,d,tab,id="none",jj=1,kill=FALSE,pause=FALSE){
  if (class(mod)=="nls"){
    Xg <- c(0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
    nX=length(Xg)
    if (DNA=="wt") {
      tinh=list(list(t=0,c=c(0,1)),list(t=1,c=0),list(t=2,c=0),list(t=10,c=0))
      cinh=list(list(t=0,c=c(0,1,5)),list(t=1,c=0),list(t=5,c=0)) } else { 
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
(tab=data.frame(id="XX",p=rep(0,1),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          KAt=0,KBt=0,
          KAc=0,KBc=0
))

(id="DEFG")
kt.A=.6
kt.B=1
kc.A=.53
kc.B=.32
KAT=2.9
KBT=106
KAC=8.5
KBC=8.2
if (id=="DEFF") {
KAC=8.4
KBC=KAC
}

mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
           +(kt.B*T/KBT/(1+T/KBT+C/KAC+t/KBt+c/KBc)))*It/2
          +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
           +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
    start=list(
        KAt=1,KBt=10,
        KAc=1,KBc=.8),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=TRUE)
print(s<-summary(mod))
s$parameters[,"Estimate"]/c(KAT,KBT,KAC,KBC)

(id="DEFG")

#0.44(0.01)  0.98(0.18)  kt.A  kt.A  2.77(0.66)  41.22(12.35)  13.49(3.38) 5.17(0.91)  0.59(0.09)
kt.A=.44
kt.B=.98
kc.A=kt.A
kc.B=kt.A
KAT=2.8
KBT=41
KAC=13.5
KBC=5.2
p=.59
if (id=="DEFF") {
  KAC=8.4
  KBC=KAC
}

mod<-nls(k~(p*(kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(1-p)*(kt.B*T/KBT/(1+T/KBT+C/KAC+t/KBt+c/KBc)))*It/2
        +(p*(kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(1-p)*(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
    start=list(
        KAt=1,KBt=10,
        KAc=1,KBc=.8),
    data=d,weights=rep(1,length(k)),model=TRUE,trace=TRUE)
tab=plotmod(mod,d,tab,id,kill=TRUE)
print(s<-summary(mod))
s$parameters[,"Estimate"]/c(KAT,KBT,KAC,KBC)



library(hwriter)
hwrite(tab,"C:/Users/radivot/case/active/papers/TK2/tc2Mono.html")

