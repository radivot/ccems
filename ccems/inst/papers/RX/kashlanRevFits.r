rm(list=ls(all=TRUE))  
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}
# read in the top models
load("/users/radivot/results/RXF3top100.RData") 
#load("/users/radivot/results/RGtXF3noHighATPtop100.RData") 
mkHTML(globalTopN)  # dump out to make sure it's what you think it is
globalTopN[[1]]
d=globalTopN[[1]]$d # grab data from first model
if (d$XF[1]<10) d=d[-1,]  # kill first point strictly to plot logs on x axis
d
lgx=log(d$XF)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(RT = rep(d$RT[1],length(fineX)), XF = fineX)
newPnts
MDS=NULL
for (i in 1:3) 
  MDS[[i]] <- simulateData(globalTopN[[i]],predict=newPnts,typeYP="m")  

if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 6, height = 6,restoreConsole = TRUE) 
par(mfcol=c(1,1),mar=c(5,4,1,1)+.1)
plot(d$XF,d$m,type="p",pch=1, xlab=expression("[ATP] ("*mu*"M)"), ylab="Mass (kDa)",log="x")
#text(100,520,"A")
for (i in 1:3)
  lines(MDS[[i]]$predict$XF,MDS[[i]]$predict$EY,type="l",lty=i)
getID <- function(x) { paste(row.names(x$params)[x$params$opt],collapse=".")}
mnms=sapply(globalTopN,getID)
#mnms[1:3]
#legend(200,520,c("Best", "Second Best", "Third Best"),lty=1:3)
legend(100,520,mnms[1:3],lty=1:3,cex=0.9,bty="n")
#plot(globalTopN[[1]]$d$EY,globalTopN[[1]]$res,type="p",pch=1, xlab="Fitted Value", ylab="Residual")
#text(200,10,"B")



