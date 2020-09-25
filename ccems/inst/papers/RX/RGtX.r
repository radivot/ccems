rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}
print(host)
  topology <- list(
      heads=c("R1X0","R2X0","R3X0"), # s-sites are already filled by dTTP 
      sites=list(                    
          a=list(                                                              # a-site       thread #
              m=c("R1X1","R1X2"),                                                 # monomer          1
              d=c("R2X1","R2X2","R2X3","R2X4"),                                          # dimer            2
              t=c("R3X1","R3X2","R3X3","R3X4","R3X5","R3X6")                            # trimer         3
          ),
          h=list( ## tails of a-site threads are heads of h-site threads       # h-site
              m=c("R1X3","R1X4"),                                                 # monomer          5
              d=c("R2X5","R2X6","R2X7","R2X8"),                                         # dimer            6
              t=c("R3X7","R3X8","R3X9","R3X10","R3X11","R3X12")                       # tetramer         7
          )
      )
  )
  dd=subset(RNR,(year==2002)&(fg==5)&(X>0)&(m>0),select=c(R,X,m,year))
dd
  names(dd)[1:2]=c("RT","XT")
#  names(dd)[1:2]=c("RT","XF")
  dd[,1]=dd[,1]/2
#  dd[1,2]=0
  #dd=dd[-12,]  # outlier that raises things to cause the downturn
#dd=dd[-14,]  # causes downturn    
#dd=dd[-c(1,12),]  # kill speculative bias at 90 and kill high at 540
#dd=dd[-c(7,8,9),]  # remove bogus first data point since 90 is built into model anyway, or it should be estimated freely
dd                   # also take out above physiological 3mM values, i.e. values at 5, 7 and 10 mM
#  
#g <- mkg(topology,free=TRUE)
g <- mkg(topology)
g

bigNms=c("localhost")
rackLen=c(12,12)
for (rack in 0:0) for (node in 0:rackLen[rack+1])  bigNms=c(bigNms,paste("compute",rack,node,sep="-"))
big=rep(4,length(bigNms))
names(big)<-bigNms
big
cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,
    "compute-0-3"=4,"compute-0-4"=4,"compute-0-5"=4,"compute-0-6"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if (host=="rnrClust") cpusPerHost=cpusPerHost[1:5]
#if (host=="rnrClust") cpusPerHost=cpusPerHost[1]
#if (host=="dck") cpusPerHost=big
if (host=="dNTP") cpusPerHost=c("localhost" = 3)
if (host=="ATP") cpusPerHost=c("localhost" = 1)
#if (host=="tk2") cpusPerHost=cpusPerHost[1]
if (host=="dck") cpusPerHost=cpusPerHost[1]
print(cpusPerHost)
tops=ems(dd,g,cpusPerHost=cpusPerHost,
    doGrids=TRUE,doSpurs=TRUE,p=-1,m1=-180,maxTotalPs=3,
#    doGrids=FALSE,doSpurs=TRUE,p=0.9,m1=90,forceM1=T,forceP=T, maxTotalPs=3,
#    doGrids=TRUE,doSpurs=FALSE,p=0.9,m1=89,forceM1=T,forceP=T, maxTotalPs=4,
    ptype="SOCK",topN=100,KIC=100,transform="none") 

