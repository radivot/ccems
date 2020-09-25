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
dd=d
dd

REMOVEM=FALSE  # nls crashes below with this set false
REMOVEM=TRUE
if (REMOVEM) {
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((T>100)&(t==10)&(c==0)&(It==1))) 
}
dim(d)

od1=subset(dd,(C==5)&(It==1)) # outliers to plot X's for reviewer
od2=subset(dd,(t==1)&(It==1)) 
od3=subset(dd,(T>100)&(t==10)&(c==0)&(It==1)) 
(od=rbind(od1,od2,od3))
# it really is t=1 that has to go, and not t=0
#d
library(hwriter)
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",p=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.Af=0,kt.Bf=0,kc.Af="kt.A",kc.Bf="kt.A",
          KAT=0,KBT=0,KAC=0,KBC="KAC",
          KAt=0, KBt=0, KAc=0, KBc=0
      ))

(id="DEFG:DE.DE.DEFG")
(id="Eq.1")
mod<-nls(k~((kt.Af*T/KAT/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(kt.Bf*T/KBT/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*It
        +((kc.Af*C/KAC/(1+T/KAT+C/KAC+t/KAt+c/KAc))
          +(kc.Bf*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)))*Ic, 
    start=list(
        kt.Af=.15, 
        kt.Bf=.5,
        kc.Af=.1,
        kc.Bf=.3,
        KAT=.46,
        KBT=106,
        KAC=8.5,
        KBC=8.2,
        KAt=1,
        KBt=10,
        KAc=1,
        KBc=.8),
    data=d,weights=rep(1,length(k)),trace=TRUE)
mod

tab=plotmodAll(mod,d,tab,id,od=od,kill=TRUE)




s<-summary(mod)$parameters[,"Estimate"]
tab
s


g=NULL
g$d=d
g$CD=FALSE
g$DNA="wt"
g$params=data.frame(initial=s,final=s,opt=TRUE,constr="none") # refit model
g=getEk(g,init=TRUE)
reportg<-function(g) {
  tmp=format(cbind(g$params[g$params$opt,"final",drop=F],g$CI),digits=2,trim=T,scientific=FALSE) 
  tmp=transform(tmp,Wald.95.CI=mapply(paste,"(",lower,", ",upper,")",sep=""))
  names(tmp)[1]<-c("Estimate")
  subset(tmp,select=c("Estimate","Wald.95.CI"))
}

g$params[,"opt"]=TRUE  
g=fitModel(g)
g$ciARD=reportg(g)
hwrite(g$ciARD,"C:/Users/radivot/case/active/papers/TK2/ciARD12.html",cellspacing=0,cellpadding=2)
gARD=g
save(gARD,file="C:/Users/radivot/case/active/papers/TK2/ARD12.RData")

nd=transform(g$d,frcA=EkA/Ek)
nd=subset(nd,select=c("T","C","t","c","It","Ic","frcA"))
nd

plotFracA(nd)

