\name{TK2}
\alias{TK2}
\docType{data}
\title{ Thymidine Kinase 2 Data}
\description{
  Human thymidine kinase 2 data.   
}
\usage{data(TK2)}
\format{
  A data frame with the following columns.
  \describe{
    \item{\code{dT}}{The total concentration of dT (nucleoside substrate).}
    \item{\code{dC}}{The total concentration of dC (nucleoside substrate).}
    \item{\code{ATP}}{The total concentration of ATP.}
    \item{\code{dTTP}}{The total concentration of dTTP.}
    \item{\code{dCTP}}{The total concentration of dCTP.}
    \item{\code{VdT}}{The velocity of the dT kinase reaction.}
    \item{\code{VdC}}{The velocity of the dC kinase reaction.}
    \item{\code{fg}}{The figure number in the original paper.}
    \item{\code{year}}{The year of the original paper.}
    \item{\code{jrnl}}{The journal of the data source: 
       \code{FEB} is FEBS Letters, 
       \code{JBC} is Journal of Biological Chemistry.
       }
    \item{\code{vol}}{The volume number of the article.}
    \item{\code{page}}{The page number of the article.}
    \item{\code{frstAut}}{The first author of the article. This is a factor.  }
    \item{\code{seq}}{Sequence of DNA. Wild type (wt) or H121N.}
    \item{\code{state}}{In Barroso papers, dimer or tetramer, i.e. separate data for these two fractions.}
    \item{\code{E}}{Enzyme concentration in uM.  }
    \item{\code{k}}{Activity in 1/sec.  }
  }
}

\details{All concentrations are in micromolar.  V is in units/mg  where
         1 unit = 1 nmol/min. Since 1000 nmoles = 29 mg, 1 mg = 1000/29 nmoles = 34.4 nmoles.
         Thus, k in 1/sec equals V/60/34.4. Enzyme concentrations can be assumed to be
         much less than the smallest total substrate concentration. For the 1991 paper, 2 ng of 
         enzyme was placed in a total reaction vessel volume of 50 uL, so the total enzyme 
         concentration was (0.002 ug/50 uL)*(1 nmol/29 ug) = 1.37e-6 nmoles/uL = 1.37 nM. 
         In the other two papers there was first a 46-fold dilution of [E] at .25 mg/mL to 5.4 ug/mL, 
         and 1 uL of this (i.e. 5.4 ng) was then placed into a 50 uL vessel. Thus, the total enzyme
         concentration in these papers was 0.0054 ug/(29 ug/nmol)/(50 uL) = 3.75e-6 nmol/uL = 3.75 nM. 
           }
\source{ The figure number, year, first author, journal, volume and page 
           number are all included in the dataframe. }
\examples{

library(ccems)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
  if (.Platform$OS.type=="windows") 
    windows(width = 5, height = 5,restoreConsole = TRUE) else 
    X11(width=5,height=5)
  par(mfcol=c(2,2),mar=c(4,4,2,1)+.1)
  for (yr in c(1991,1999)) {
  d=subset(TK2,(dT>0)&(year==yr),select=c(dT,VdT,frstAut,year))
  d=transform(d,V=VdT,S=dT)
  plot(d$S,d$V,xlab="Total [dT]",log="xy", ylab="V ", 
        main=paste(d[1,"frstAut"],d[1,"year"]))
  d=subset(TK2,(dC>0)&(year==yr),select=c(dC,VdC,frstAut,year))
  d=transform(d,V=VdC,S=dC)
  plot(d$S,d$V,xlab="Total [dC]",log="xy", ylab="V ", 
        main=paste(d[1,"frstAut"],d[1,"year"]))
  }
 }
\keyword{datasets}
