\name{TK2BBA77}
\alias{TK2BBA77}
\docType{data}
\title{ Thymidine Kinase 2 Data that underlies Biochimica et Biophysica Acta. 478 (1977) 364-375 }
\description{
 This is the healthy donor lymphocyte thymidine kinase 2 data that underlies
 Biochimica et Biophysica Acta. 478 (1977) 364-375 (thus the name TK2BBA77).  
}
\usage{data(TK2BBA77)}
\format{
  A data frame with the following columns.
  \describe{
    \item{\code{S}}{The total concentration of dT (thymidine).}
    \item{\code{V}}{Velocity (arbitrary units since concentrations are unknown).}
    \item{\code{donor}}{Healthy blood donor number. Lymphocytes of donors 1-4 were PHA stimulated, the rest
                    were not and thus represent G0 cells. }
  }
}

\details{All concentrations are in micromolar. Enzyme concentrations can be assumed to be
         much less than the smallest total substrate concentration.   }
\source{ Munch-Petersen and Tyrsted, Biochimica et Biophysica Acta. 478 (1977) 364-375. }
\examples{

library(ccems)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
  windows(width = 9, height = 6,restoreConsole = TRUE,ypos=0) else X11(width=9,height=6)
par(mfrow=c(2,4),mar=c(4,4,2,1)+.1)
for (j in 1:8){
  d=subset(TK2BBA77,donor==j)
  print(d)
  try(MMmod<-nls(V~Va*S/(S50a+S)+Vb*S/(S50b+S), 
      start=list(Va=.15,Vb=.15,S50a=1,S50b=15),data=d,weights=1/V^2,algorithm="port",lower=0))
  print(summary(MMmod))
  plot(d$S,d$V,xlab="Total [dT]",log="xy", ylab="velocity (arb. units) ", 
      main=paste("Donor #",d[1,"donor"],sep="") )
  lgx=log(d$S)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  lines(fineX,predict(MMmod,list(S=fineX)),col="black",lwd=1)
  mtext(paste("S50a = ",format(MMmod$m$getPars()["S50a"],digits=3),sep=""),line=-2.2,side=1,font=1,cex=0.7)
  mtext(paste("S50b = ",format(MMmod$m$getPars()["S50b"],digits=3),sep=""),line=-1.2,side=1,font=1,cex=0.7)
}

 }
\keyword{datasets}
