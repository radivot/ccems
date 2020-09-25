rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK2/data/TK2.rda",sep=""))
}
TK2=transform(TK2,ET=E/2)
TK2
d=subset(TK2,(dCTP==0)&(dTTP==0)&(ATP>=2000)&(year==2003)&(seq=="wt"),select=c(ET,dT,dC,k,VdT,VdC,seq,frstAut,year))
names(d)[2:3]=c("TF","CF")
dt=d[!is.na(d$VdT),]
dc=d[!is.na(d$VdC),]
dt
dc

topology <- list(  
    heads=c("E1S0","E2S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site 
            m=c("E1S1"),
            d= c("E2S1","E2S2")  
        ) # d for dimer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=TRUE,free=TRUE)
#g <-mkg(topology, activity=TRUE,TCC=TRUE)
g
gg=mkGrids(g, maxTotalPs=3,m1=-30)  # TK2 is 30 kDa, though irrelevant here
gg$chunk


tops=ems(d,g,
    doGrids=TRUE,doSpurs=FALSE,p=-1,m1=-29,forceM1=F,forceP=F, maxTotalPs=2,
#    doGrids=FALSE,doSpurs=TRUE,p=0.9,m1=90,forceM1=T,forceP=T, maxTotalPs=3,
#    doGrids=TRUE,doSpurs=FALSE,p=0.9,m1=89,forceM1=T,forceP=T, maxTotalPs=4,
    topN=100,KIC=1) 


#gg=mkSpurs(g)
#gg$chunk
#
#mkGrids(g)
#

cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,
    "compute-0-3"=4,"compute-0-4"=4,"compute-0-5"=4,"compute-0-6"=4) # for tk2
if (host=="rnrClust") cpusPerHost=cpusPerHost[1:5]
if (host=="ATP") cpusPerHost=c("localhost" = 1)
print(cpusPerHost)
tops=ems(d,g,cpusPerHost=cpusPerHost,
    doGrids=FALSE,doSpurs=TRUE,p=-1,m1=-29,forceM1=F,forceP=F, maxTotalPs=2,
#    doGrids=FALSE,doSpurs=TRUE,p=0.9,m1=90,forceM1=T,forceP=T, maxTotalPs=3,
#    doGrids=TRUE,doSpurs=FALSE,p=0.9,m1=89,forceM1=T,forceP=T, maxTotalPs=4,
    ptype="SOCK",topN=100,KIC=1,transform="none") 

### Warning: the next line clears all existing figures!!
#if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
#if (.Platform$OS.type=="windows") 
#  windows(width = 3, height = 6,restoreConsole = TRUE,ypos=0) else X11(width=2,height=4)
#par(mfcol=c(2,1),mar=c(4,4,2,1)+.1)
