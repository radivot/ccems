# in Discussion
B=c(1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975)
y=choose(10,1:10)
sum(B*y)

n=1:6
(1+(1-1/n)^0.5)/(1-(1-1/n)^0.5)

#####

rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/ccems/ccems/inst/papers/TK2/data/TK2.rda",sep=""))
}

topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)


KS=c("E1S0_S","E1S1_S")
kS=c("kE1S1","kE1S2")

options(digits=2)

fineX=seq(1,25,by=1)
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
a=5
b=0.5
testCases=list(
    list(Kic=a*c( 1,1), kic=b*c(1.5,.67),hinit=2,kminit=.67*b),
    list(Kic=a*c( 1,1), kic=b*c(1.5,.67)[2:1],hinit=1.1,kminit=1.5*b),
    list(Kic=a*c( 1,1), kic=b*c(3,.33),hinit=1.3,kminit=1*b),
    list(Kic=a*c( 1,1), kic=b*c(3,.33)[2:1],hinit=1.3,kminit=3*b),
    list(Kic=a*c(1.5,.67)[2:1], kic=b*c(1,1),hinit=1,kminit=1*b),
    list(Kic=a*c(1.5,.67), kic=b*c(1,1),hinit=1,kminit=1*b),
    list(Kic=a*c(3.3,.33)[2:1], kic=b*c(1,1),hinit=1,kminit=1*b),
    list(Kic=a*c(3.3,.33), kic=b*c(1,1),hinit=1,kminit=1*b)
#    list(Kic=c(16,.3)[2:1], kic=b*c(1,1),hinit=0.6,kminit=1*b),
#    list(Kic=c(16,.3), kic=b*c(1,1),hinit=1,kminit=1*b)
)

Kmapping = mkKd2Kj(g)
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 10, height = 5,restoreConsole = TRUE,ypos=-25) else X11(width=12.5,height=5)
par(mfcol=c(2,4),mar=c(4,4,2,1)+.1)

for (i in 1:length(testCases)) {
  Kic=testCases[[i]]$Kic  
  kic=testCases[[i]]$kic  
  hinit=testCases[[i]]$hinit  
  kminit=b*testCases[[i]]$kminit  
  names(Kic)<-KS
  names(kic)<-kS
  Kstrn="         K = "
  for (j in KS) 
    Kstrn=paste(Kstrn,sprintf("%3.1f",Kic[j]),sep=ifelse(j!="E1S0_S",", ",""))
  Kstrn 
  kstrn="         k = "
  for (j in kS) 
    kstrn=paste(kstrn,sprintf("%3.1f",kic[j]),sep=ifelse(j!="kE1S1",", ",""))
  kstrn
  mdl = mkModel(g,"full",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  names(df)[3]<-"k"
  plot(df$ST,df$k,xlab=expression(group("[",S[T],"]")),ylab="k")
  mtext(Kstrn,line=-1.2,side=1,font=1,cex=0.7,adj=.45)
  mtext(kstrn,line=-2.2,side=1,font=1,cex=0.7,adj=.45)
  print(Kstrn)
  print(kstrn)
  e<-try(hillda<-nls(k~kmax*(ST/S50)^h/(1+(ST/S50)^h),df,algorithm="port",lower=c(.1,1,.6),
          upper=c(5,20,3),start=list(kmax=kminit,S50=5,h=hinit)))
  if (class(e)=="try-error") 
  {print("******************************  FIT FAILED ****************")} else
  {
    print(hillda)
    lines(fineX,predict(hillda,list(ST=fineX)),col="black",lwd=1)
    mtext(bquote(k[max] == .(format(hillda$m$getPars()["kmax"],digits=2))),line=-6.2,side=1,font=1,cex=0.7,adj=.9)
    mtext(bquote(S[50] == .(format(hillda$m$getPars()["S50"],digits=2))),line=-5.2,side=1,font=1,cex=0.7,adj=.9)
# mtext(paste("kmax = ",format(hillda$m$getPars()["kmax"],digits=2),sep=""),line=-6.2,side=1,font=1,cex=0.7,adj=.9)
# mtext(paste("S50 = ",format(hillda$m$getPars()["S50"],digits=2),sep=""),line=-5.2,side=1,font=1,cex=0.7,adj=.9)
    mtext(paste("h = ",format(hillda$m$getPars()["h"],digits=2),sep=""),line=-4.2,side=1,font=1,cex=0.7,adj=.9)
  }
  
  
}

#i=7  # check to see if correlated residuals are still obvious in a hillplot. Answer: yes 
#Kic=testCases[[i]]$Kic  
#kic=testCases[[i]]$kic  
#hinit=testCases[[i]]$hinit  
#kminit=b*testCases[[i]]$kminit  
#names(Kic)<-KS
#names(kic)<-kS
#Kstrn="         K = "
#for (j in KS) 
#  Kstrn=paste(Kstrn,sprintf("%3.1f",Kic[j]),sep=ifelse(j!="E1S0_S",", ",""))
#Kstrn 
#kstrn="         k = "
#for (j in kS) 
#  kstrn=paste(kstrn,sprintf("%3.1f",kic[j]),sep=ifelse(j!="kE1S1",", ",""))
#kstrn
#mdl = mkModel(g,"full",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
#df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
#names(df)[3]<-"k"
#hillda<-nls(k~kmax*(ST/S50)^h/(1+(ST/S50)^h),df,algorithm="port",lower=c(.1,1,.6),
#        upper=c(5,20,3),start=list(kmax=kminit,S50=5,h=hinit))
#Ek=predict(hillda,list(ST=df$ST))
#s<-summary(hillda)
#kmax=s$parameters["kmax","Estimate"]
#windows(width = 5, height = 5,restoreConsole = TRUE,ypos=-25) 
#par(mfcol=c(1,1),mar=c(4,4,2,1)+.1)
#
#plot(df$ST,df$k/(kmax-df$k),log="xy",xlab=expression(group("[",S[T],"]")),ylab="k")
#lines(df$ST,Ek/(kmax-Ek))
#



