rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
TK2

DNA="H121N"
DNA="wt"
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[1:2]=c("T","C")
dtm=subset(d,(It==1)&(C==0))
dcm=subset(d,(Ic==1)&(T==0))
dm=rbind(dtm,dcm) # m for marginals
#if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
#if (.Platform$OS.type=="windows") 
#  windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
#par(mfcol=c(1,2),mar=c(4.2,4,1,1)+.1)
#with(dtm,plot(T,k,log="x", pch=as.numeric(as.factor(year)),col=as.numeric(as.factor(year)),
#        ylab="dT Kinase activity (1/sec)",xlab="[dT] uM" ))
#legend("topleft",legend=c("1991","1999","2003"), pch=1:3,col=1:3,bty="n")
#with(dcm,plot(C,k,log="x", pch=as.numeric(as.factor(year)),col=as.numeric(as.factor(year)),
#        ylab="dC Kinase activity (1/sec)",xlab="[dC] uM" ))

d1=subset(dm,year==1991)
d2=subset(dm,year==1999)
d3=subset(dm,(year==2003)&(frstAut=="Wang"))
years=c(1991,1999,2003)

plotmods<-function(mods,ds,ids,picks,kill=FALSE){
  Xg <- c(.075,.1,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
  pd0=rbind(data.frame(T=Xg,C=0,It=1,Ic=0),
      data.frame(T=0,C=Xg,It=0,Ic=1) )
  if (kill) if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
  if (.Platform$OS.type=="windows") 
    windows(width = 8, height = 4,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
  par(mfcol=c(1,2),mar=c(4.2,4,0,1)+.1,oma=c(0,0,0.5,0))
  plot(Xg,type="n",ylim=range(c(0,ds$k)),xlim=range(Xg),ylab="dT Kinase activity (1/sec)",xlab="[dT] uM",log="x")
  for (j in 1:3) 
  {
      pd=data.frame(pd0,Ek=predict(mods[[picks[j]]],pd0))
      d=mods[[picks[j]]]$d
#      d=subset(ds,year==years[j])
      with(d[d$It>0,],points(T,k,col=j,pch=j))
      with(pd[pd$It>0,],lines(T,Ek,col=j))
      print(s<-summary(mods[[picks[j]]]))
  }
  legend("topleft",legend=ids[picks], lty=1,col=1:3,bty="n")
  plot(Xg,type="n",ylim=range(c(0,ds$k)),xlim=range(Xg),ylab="dC Kinase Activity (1/sec)",xlab="[dC] uM",log="x")
  for (j in 1:3) 
    if (class(mods[[picks[j]]])=="nls"){
      pd=data.frame(pd0,Ek=predict(mods[[picks[j]]],pd0))
      d=mods[[picks[j]]]$d
#      d=subset(ds,year==years[j])
      with(d[d$Ic>0,],points(C,k,col=j,pch=j))
      with(pd[pd$Ic>0,],lines(C,Ek,col=j))
    }
  legend("topleft",legend=c("1991 Munch-Petersen","1999 Wang","2003 Wang"), pch=1:3,col=1:3,bty="n")
}

mkTab<-function(mods,ds,ids,jj=1){
options(stringsAsFactors = FALSE)
(tab=data.frame(year=rep(0,length(mods)),id="XX",p=0,SSE=0,AIC=0,kt.T=0,kt.TT="kt.T",kc.C="kt.T",kc.CC="kc.C",KT=0,KTT=0,KC=0,KCC="KC"))
  Xg <- c(.075,.1,.2,.35,0.5,0.75,1,1.5,2,3.5,5,7.5,10,15,20,35,50,100,150,200)
  pd0=rbind(data.frame(T=Xg,C=0,It=1,Ic=0),
      data.frame(T=0,C=Xg,It=0,Ic=1) )
  for (j in 1:length(mods)) 
  {
    if (class(mods[[j]])=="nls"){
      pd=data.frame(pd0,Ek=predict(mods[[j]],pd0))
      d=mods[[j]]$d
#      d=subset(ds,year==mods[[j]]$year)
      print(s<-summary(mods[[j]]))
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
      tmp=c(year=mods[[j]]$year,id=ids[j],p=p,SSE=SSE,AIC=AICc,vals)
      if (length(which(tab$id!="XX"))==0) jj=1 else jj=max(which(tab$id!="XX"))+1
      tab[jj,names(tmp)]=tmp
    } else {tab[jj,1]=ids[j]; tab[jj,3]="failed";tab[jj,4:dim(tab)[2]]="X";}
  }
  tab
}

#kc(EC) = kc(ECC) = kc(ECT)  =.21/sec, kt(ET) = kt(ETC) =.3/sec , kt(ETT) = .4/sec, 
#KE_C = KEC_C= 8.5 µM, KE_T = 5 µM and KT_T = 46 µM 

mods=NULL
ids=NULL
(ids[1]="DEFF:DE.DD")
mods[[1]]<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=.3,KTT=16,KC=5),
    data=d1,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[1]]$year=1991
mods[[1]]$d=d1

(ids[2]="DDFF:DE.DD")
mods[[2]]<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.2,kc.C=.2,KT=.3,KTT=23,KC=43),
    data=d1,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[2]]$year=1991
mods[[2]]$d=d1

(ids[3]="DDDD:DE.DD")
mods[[3]]<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,KT=.6,KTT=23,KC=43),
    data=d1,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[3]]$year=1991
mods[[3]]$d=d1


(ids[4]="DEFF:DE.DD")
mods[[4]]<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.2,kt.TT=.35,kc.C=.4,KT=0.86,KTT=8.5,KC=20),
    data=d2,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[4]]$year=1999
mods[[4]]$d=d2

(ids[5]="DEEE:DE.DD")
mods[[5]]<-
    nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kt.TT*2*C/KC+2*kt.TT*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.2,kt.TT=.35,KT=0.86,KTT=8.5,KC=20),
        data=d2,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[5]]$year=1999
mods[[5]]$d=d2

(ids[6]="DDDD:DE.DD")
mods[[6]]<-
    nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
            (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
        start=list(kt.T=.35,KT=0.86,KTT=8.5,KC=20),
        data=d2,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[6]]$year=1999
mods[[6]]$d=d2

(ids[7]="DEFF:DE.DD")
mods[[7]]<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=10,KTT=100,KC=10),
    data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[7]]$year=2003
mods[[7]]$d=d3

(ids[8]="DDFF:DE.DD")
mods[[8]]<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kc.C*2*C/KC+2*kc.C*C^2/(KC*KC))*Ic)/
   (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.8,kc.C=.6,KT=10,KTT=100,KC=10),
    data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[8]]$year=2003
mods[[8]]$d=d3

(ids[9]="DDDD:DE.DD")
mods[[9]]<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC+2*kt.T*C^2/(KC*KC))*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC+C^2/(KC*KC))/2, 
    start=list(kt.T=.52,KT=8,KTT=7,KC=16),
#    data=d3,weights=rep(1,length(k)),trace=TRUE)
data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[9]]$year=2003
mods[[9]]$d=d3

(ids[10]="DEFX:DE.DI")
mods[[10]]<-nls(k~((kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))*It+(kc.C*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.4,kt.TT=.8,kc.C=.6,KT=10,KTT=100,KC=10),
    data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[10]]$year=2003
mods[[10]]$d=d3

(ids[11]="DDDX:DE.DI")
mods[[11]]<-nls(k~((kt.T*2*T/KT+2*kt.T*T^2/(KT*KTT))*It+(kt.T*2*C/KC)*Ic)/
        (1+2*T/KT+T^2/(KT*KTT)+2*C/KC)/2, 
    start=list(kt.T=.4,KT=10,KTT=100,KC=10),
    data=d3,weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
mods[[11]]$year=2003
mods[[11]]$d=d3

tab=mkTab(mods,dm,ids)
tab
tab[tab[,"kc.C"]=="kt.T","kc.CC"]="kt.T"
tab[5,c("kc.C","kc.CC")]="kt.TT"
tab[10:11,"KCC"]="Inf"
names(tab)[2]<-"Model"
tab

library(hwriter)
hwrite(tab,"C:/Users/radivot/case/active/papers/TK2/tab4.html",cellspacing=0,cellpadding=2)

aics=as.numeric(tab$AIC)
mins=tapply(aics,tab$year,min)
picks=NULL
for (i in 1:length(mins)) picks[i]=which(aics==mins[i])
picks
plotmods(mods,dm,ids,picks,kill=T)


#nls(k~(kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))/
#        (1+2*T/KT+T^2/(KT*KTT)), 
#    start=list(kt.T=.1,kt.TT=.25,KT=0.5,KTT=15),
##    start=list(kt.T=.09,kt.TT=.17,KT=0.86,KTT=8.5),
##    data=d2[d2$It==1,],weights=rep(1,length(k)),trace=TRUE)
#    data=d2[d2$It==1,],weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)
#
#nls(k~(kt.T*2*T/KT+2*kt.TT*T^2/(KT*KTT))/
#        (1+2*T/KT+T^2/(KT*KTT)), 
#    start=list(kt.T=.1,kt.TT=.25,KT=0.5,KTT=15),
##    data=dd,weights=rep(1,length(k)),trace=TRUE)
#    data=d2[d2$It==1,],weights=rep(1,length(k)),algorithm="port",lower=0,trace=TRUE)

