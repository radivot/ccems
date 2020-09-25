rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,It=as.numeric(!is.na(TK2$VdT)),Ic=as.numeric(!is.na(TK2$VdC)))
DNA="wt"
d=subset(TK2,(year==2003)&(frstAut=="Wang")&(ATP>=2000)&(seq==DNA))#,select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[c(1,2,4,5)]=c("T","C","t","c")
dd=d
dd
d=subset(d,!((C==5)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((t==1)&(It==1))) # kill outlier profile which accounts for 75% of SSE 
d=subset(d,!((T>100)&(t==10)&(c==0)&(It==1))) 
dim(d)
od1=subset(dd,(C==5)&(It==1)) # outliers to plot X's for reviewer
od2=subset(dd,(t==1)&(It==1)) 
od3=subset(dd,(T>100)&(t==10)&(c==0)&(It==1)) 
(od=rbind(od1,od2,od3))

library(hwriter)
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",p=rep(0,6),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          h=0,kt.A=0,kc.B=0,
          KAT=0,KBT=0,KAC=0,KBC=0,
          KAt=0, KBt=0, KAc=0, KBc=0
      ))

(id="Eq. 2")

mod<-nls(k~((kt.A*(T/KAT)^h/(1+(T/KAT)^h+t/KAt+c/KAc)/(1+C/KAC)))*It
          +((kc.B*C/KBC/(1+T/KBT+C/KBC+t/KBt+c/KBc)/(1+t/KBt)))*Ic, 
    start=strt<-list(
        h=0.5,
        kt.A=.25, 
        kc.B=.8,
        KAT=13,
        KBT=4.9,
        KAC=40,
        KBC=11,
        KAt=2,
        KBt=2.5,
        KAc=.79,
        KBc=.87),
    data=d,weights=rep(1,length(k)),trace=TRUE)
tab=plotmodAll(mod,d,tab,id,od=od,kill=TRUE)

mod
ci95=confint(mod)
ci95
s1=summary(mod)
str(mod)
#cbind(strt,s1$coefficients[,"Estimate",drop=F],ci95)
data.frame(start=t(t(strt)),s1$coefficients[,"Estimate",drop=F],ci95)
#hwrite(df,"C:/Users/radivot/case/active/papers/TK2/eq2.html",cellspacing=0,cellpadding=2)



s<-summary(mod)$parameters[,"Estimate"]
tab
s


g=NULL
g$d=d
g$CD=TRUE
g$DNA="wt"
g$params=data.frame(initial=t(t(strt)),final=s,opt=TRUE,constr="none") # refit model
g

g=getEk(g,init=TRUE)


g$params[,"opt"]=TRUE  
g=fitModel(g)
reportg<-function(g) {
  tmp=format(cbind(g$params[g$params$opt,c("initial","final"),drop=F],g$CI),digits=2,trim=T,scientific=FALSE) 
  tmp=transform(tmp,Wald.95.CI=mapply(paste,"(",lower,", ",upper,")",sep=""))
  names(tmp)[2]<-c("Estimate")
  subset(tmp,select=c("initial","Estimate","Wald.95.CI"))
}
g$ciCD=reportg(g)
hwrite(g$ciCD,"C:/Users/radivot/case/active/papers/TK2/ciCD12.html",cellspacing=0,cellpadding=2)
gCD=g
save(gCD,file="C:/Users/radivot/case/active/papers/TK2/CD12.RData")


