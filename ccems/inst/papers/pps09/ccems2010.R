rm(list=ls()) #this is pps09.R with minimal local computing of 1-param models
# library(tools)
# showNonASCIIfile("~/case/active/jarek/Rjarek/p53/ccems2010.R")
# use https://pages.cs.wisc.edu/~markm/ascii.html
library(ccems) 
topology <- list(
  heads=c('R1X0', 'R2X2', 'R4X4', 'R6X6'),
  sites=list(
    a=list(
      m=c('R1X1'), 
      d=c('R2X3', 'R2X4'),
      t=c('R4X5', 'R4X6', 'R4X7', 'R4X8'),
      h=c('R6X7', 'R6X8', 'R6X9', 'R6X10', 'R6X11', 'R6X12') 
    ), 
    h=list(
      m=c('R1X2'), 
      d=c('R2X5', 'R2X6'), 
      t=c('R4X9', 'R4X10', 'R4X11', 'R4X12'),
      h=c('R6X13', 'R6X14', 'R6X15', 'R6X16', 'R6X17', 'R6X18') 
    )
  )
)  
g=mkg(topology)
dd=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R, X, m, year))
names(dd)[1:2]=c('RT', 'XT') 
cpus=c('localhost'=4)
top10=ems(dd,g,cpusPerHost=cpus,maxTotalPs=1,ptype='SOCK',
          KIC=100, topN=10, transform= 'none')
# Time difference of 4.533714 mins
# Fitted =  13 , out of a total of  13 

