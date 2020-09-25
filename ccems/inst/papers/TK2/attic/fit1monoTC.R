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
library(hwriter)
source(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/common.R",sep=""))

# ONE MONOMER Model
options(stringsAsFactors = FALSE)
(tab=data.frame(id="XX",p=rep(0,1),SSE=0,AIC=0, #sig=0,AICc=0,AIC=0,AICchk=0,
          kt.T=0,kc.C="kt",KT=0,KC="KT"))

mod<-nls(k~(kt.T*T/KT/(1+T/KT+C/KC))*It
          +(kc.C*C/KC/(1+T/KT+C/KC))*Ic, 
    start=list(
        kt.T=.6,
        kc.C=.6,
        KT=5,
        KC=7),
    data=d,weights=rep(1,length(k)),trace=TRUE)
(mod$id="DE:DE")
mod$data="2003 Wang"
mod$d=d
tab=plotmod(mod,tab,kill=T)
summary(mod)  # note that the K CI overlap and can thus be pooled

mod<-nls(k~(kt.T*T/KT/(1+T/KT+C/KT))*It
          +(kc.C*C/KT/(1+T/KT+C/KT))*Ic, 
    start=list(
        kt.T=.6,
        kc.C=.6,
        KT=7),
    data=d,weights=rep(1,length(k)),trace=TRUE)
(mod$id="DE:DE")
mod$data="2003 Wang"
mod$d=d
tab=plotmod(mod,tab)
summary(mod) # the k sd are 10 fold less than delta mu, so pooling should fail, as it does next

mod<-nls(k~(kt.T*T/KC/(1+T/KC+C/KC))*It
          +(kt.T*C/KC/(1+T/KC+C/KC))*Ic, 
    start=list(
        kt.T=.6,
        KC=7),
    data=d,weights=rep(1,length(k)),trace=TRUE)
(mod$id="DD:DD")
mod$data="2003 Wang"
mod$d=d
tab=plotmod(mod,tab)
tab1=tab
tab1
hwrite(tab1,"C:/Users/radivot/case/active/papers/TK2/oneMonomer.html")

