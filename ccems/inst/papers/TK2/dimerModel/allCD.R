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
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((T>100)&(t==10)&(c==0)&(It==1))) 
dim(d)

# it really is t=1 that has to go, and not t=0
#d
library(hwriter)
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",p=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
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

#(id="DEFG:DE.DE.DEFF")
#mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KAc)))*It/2
#        +((kc.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KAc)))*Ic/2, 
#    start=list(
#        kt.A=.2, 
#        kt.B=1,
#        kc.A=.2,
#        kc.B=.6,
#        KAT=.46,
#        KBT=106,
#        KAC=8.5,
#        KBC=8.2,
#        KAt=1,
#        KBt=10,
#        KAc=1),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmodAll(mod,d,tab,id,kill=TRUE)
#
#
#(id="DEDG:DE.DE.DEFG")
#mod<-nls(k~((kt.A*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kt.B*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It/2
#        +((kt.A*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
#          +(kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic/2, 
#    start=list(
#        kt.A=.2, 
#        kt.B=1,
##        kc.A=.2,
#        kc.B=.6,
#        KAT=.46,
#        KBT=106,
#        KAC=8.5,
#        KBC=8.2,
#        KAt=1,
#        KBt=10,
#        KAc=1,
#        KBc=.8),
#    data=d,weights=rep(1,length(k)),trace=TRUE)
#tab=plotmodAll(mod,d,tab,id,kill=TRUE)




#c(0.233, 0.978, 0.305, 0.613, 1.112, 17.95, 30.84,  6.200    0.003    6.916   13.851    0.624
s<-summary(mod)$parameters[,"Estimate"]
ardP=as.list(s)
cd24=ARD12(ardP)
#library(hwriter)
#cd24df=as.data.frame(lapply(cd24,format,digits=2))
#hwrite(cd24df,"C:/Users/radivot/case/active/papers/TK2/ARDall.html",cellspacing=0,cellpadding=2)
#cd24df

g=NULL
g$d=d
g$CD=FALSE
g$CD=TRUE
g$DNA="wt"
if (g$CD) {  g$params=data.frame(initial=unlist(cd24),final=unlist(cd24),opt=TRUE,constr="none")
}  else {  g$params=data.frame(initial=s,final=s,opt=TRUE,constr="none") }# else refit the ARD model

g=getEk(g,init=TRUE)

reportg<-function(g) {
 tmp=format(cbind(g$params[g$params$opt,"final",drop=F],g$CI),digits=2,trim=T,scientific=FALSE) 
 tmp=transform(tmp,Wald.95.CI=mapply(paste,"(",lower,", ",upper,")",sep=""))
 names(tmp)[1]<-c("Estimate")
 subset(tmp,select=c("Estimate","Wald.95.CI"))
}

if (!g$CD) {
  g$params[,"opt"]=TRUE  
  g=fitModel(g)
  g$ciARD=reportg(g)
  hwrite(g$ciARD,"C:/Users/radivot/case/active/papers/TK2/ciARD12.html",cellspacing=0,cellpadding=2)
  gARD=g
  save(gARD,file="C:/Users/radivot/case/active/papers/TK2/ARD12.RData")
  
} else {
  
#gg$params[,"opt"]=TRUE  
# this block tells us that g24 is the starting point
# a similar block proved that no single parameter changes would improve the fit
#mods=NULL
#mod=NULL
#for (i in 24) {
#  g=gg
#  g$params[i,"opt"]=FALSE
#  e=try(g<-fitModel(g))
#  if (class(e)=="try-error") mod="failed" else {
#    mod$params=g$params
#    mod$SSE=g$SSE
#    mod$AIC=g$AIC
#  }
#  mods[[i]]=mod
#}
#mods

g$params[,"opt"]=TRUE  
g$params[24,"opt"]=FALSE
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=fitModel(g)
g24=g

options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",p=0,SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,
          kt.TT=0,
          kt.TC=0,
          kt.tT="kt.TT",
          kt.cT=0,
          kc.C="kt.TT",
          kc.CC="kt.TT",
          kc.TC=rep("kt.T",5),
          kc.tC="kt.T",
          kc.cC="kt.TC  ",
          KT=0,
          KTT=0,
          KtT=0,
          KcT="KtT",
          KC=0,
          KCC=0,
          KTC="KC",
          KtC=c(0,0,0,0,0),
          KcC="KCC",
          Kt=0,
          Kc=0,
          Ktt=0,
          Ktc=Inf,
          Kcc=Inf))
g24$id="DEFGK,LMNRS:DEFG.DEFGK.DEFGk"
tab=plotAllg(g24,tab)

g=g24
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kc.C","kc.CC","kc.TC","kc.tC","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("Kcc"),"constr"]="Ktc"
g=fitModel(g)
g
g4678924=g
g4678924$id="DEFEK,EDDDS:DEFG.DEFGK.DEFGG"
tab=plotAllg(g4678924,tab)

g=g4678924
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("Kcc"),"constr"]="Ktc"
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=fitModel(g)
g
g467891024=g

g467891024$id="DEFEK,EDDDF:DEFG.DEFGK.DEFGG"
tab=plotAllg(g467891024,tab)

gCD=g467891024
save(gCD,file="C:/Users/radivot/case/active/papers/TK2/CD17.RData")



g=g467891024
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kt.cT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params["kt.cT",c("initial","final")]=0
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("Kcc"),"constr"]="Ktc"
g=fitModel(g)
g
g4567891024=g

g4567891024$id="DEFE0,EDDDF:DEFG.DEFGK.DEFGG"
tab=plotAllg(g4567891024,tab)

g=g4567891024
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g     # all else is stable, so set these Ktc and Kcc to Inf
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kt.cT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","Ktc","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params["kt.cT",c("initial","final")]=0
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("Ktc","Kcc"),c("initial","final")]=Inf
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g
g45678910II=g

g45678910II$id="DEFE0,EDDDF:DEFG.DEFGK.DEFII"
tab=plotAllg(g45678910II,tab)


g=g45678910II
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kt.cT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","KTC","Ktc","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params["kt.cT",c("initial","final")]=0
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("KTC"),"constr"]="KC"
g$params[c("Ktc","Kcc"),c("initial","final")]=Inf
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g
g4567891017II=g

g4567891017II$id="DEFE0,EDDDF:DEFG.DEDGK.DEFII"
tab=plotAllg(g4567891017II,tab)

g=g4567891017II
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kt.cT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","KTC","KcC","Ktc","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params["kt.cT",c("initial","final")]=0
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("KTC"),"constr"]="KC"
g$params[c("KcC"),"constr"]="KCC"
g$params[c("Ktc","Kcc"),c("initial","final")]=Inf
g=fitModel(g)
g
g456789101719II=g

g456789101719II$id="DEFE0,EDDDF:DEFG.DEDGE.DEFII"
tab=plotAllg(g456789101719II,tab)
# 1     2     3     4     5       6    7     8     9    10      11   12     13   14    15   16    17    18   19    20   21   22   23   24
#kt.T kt.TT kt.TC kt.tT kt.cT   kc.C kc.CC kc.TC kc.tC kc.cC    KT   KTT   KtT   KcT   KC   KCC   KTC   KtC  KcC   Kt   Kc   Ktt  Ktc  Kcc
#0.64  0.71  0.39  0.67  0.26   0.72  0.52  0.53  0.55  0.51   6.96 37.79  13.87 7.54 13.04 82.59 10.52 6.66 37.24 0.44 1.68 8.41 0.96 7.24
g=g456789101719II
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g$params["opt"]=TRUE
g$params[c("kt.tT","kt.cT","kc.C","kc.CC","kc.TC","kc.tC","kc.cC","KcT","KTC","KcC","Ktc","Kcc"),"opt"]=FALSE
g$params["constr"]="none"
g$params[c("kt.tT","kc.C"),"constr"]="kt.TT"
g$params["kt.cT",c("initial","final")]=0
g$params[c("kc.CC","kc.TC","kc.tC"),"constr"]="kt.T"
g$params[c("kc.cC"),"constr"]="kt.TC"
g$params[c("KTC"),"constr"]="KC"
g$params[c("KcT"),"constr"]="KtT"
g$params[c("KcC"),"constr"]="KCC"
g$params[c("Ktc","Kcc"),c("initial","final")]=Inf
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g$params["initial"]=g$params["final"]
g=getEk(g,init=TRUE)
g=fitModel(g)
g
g45678910141719II=g

g45678910141719II$id="DEFE0,EDDDF:DEFF.DEDGE.DEFII"
tab=plotAllg(g45678910141719II,tab)

reportg<-function(g) {
  tmp=format(cbind(g$params[g$params$opt,"final",drop=F],g$CI),digits=2,trim=T,scientific=FALSE) 
  tmp=transform(tmp,Wald.95.CI=mapply(paste,"(",lower,", ",upper,")",sep=""))
  names(tmp)[1]<-c("Estimate")
  subset(tmp,select=c("Estimate","Wald.95.CI"))
}
g45678910141719II$ciCD=reportg(g45678910141719II)
hwrite(g45678910141719II$ciCD,"C:/Users/radivot/case/active/papers/TK2/ciCD12.html",cellspacing=0,cellpadding=2)
gCD=g45678910141719II
save(gCD,file="C:/Users/radivot/case/active/papers/TK2/CD12.RData")

#cd24df=as.data.frame(lapply(as.data.frame(t(g$params[,"final",drop=F])),format,digits=2))
#hwrite(cd24df,"C:/Users/radivot/case/active/papers/TK2/CD24.html",cellspacing=0,cellpadding=2)
tmp=tab
tab=tmp
#tab$id<-gsub(",",",\n",tab$id,fixed=T)
tab$id<-gsub(".",".\n",tab$id,fixed=T)
tab$id<-gsub(":",":\n",tab$id,fixed=T)
names(tab)[1]<-"Model"
tab
hwrite(t(tab),"C:/Users/radivot/case/active/papers/TK2/allCDs.html",cellspacing=0,cellpadding=2)
}  # end else on g$CD = TRUE

