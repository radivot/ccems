setwd("/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/data/barroso") 
d=read.csv("barroso03dimer.csv",header=T)
#Vm=c(573,537,277,264)
nd=rbind(data.frame(T=d[1:22,1],C=0,Vm=573,y=1/d[1:22,2]),
         data.frame(T=0,C=d[23:45,1],Vm=277,y=1/d[23:45,2]))
# y=(Vm-V)/V => Vy+V=Vm => V=Vm/(y+1)
nd = transform(nd,V=Vm/(y+1))
#nd
ST=c(.05,.1,.25,.5,1,1.5,2,3,4,5,7.5,10,12.5,15,20,30,40,50,75,100,150,200)
SC=c(.05,.1,.25,.5,1,1.5,2,3,4,5,7.5,10,12.5,15,20,30,40,50,75,100,150,200,300)
dV=read.csv("barroso03Vdimer.csv",header=T) # this is crap relative to log scale data
#dV
dd=data.frame(T=c(ST,rep(0,length(SC))),C=c(rep(0,length(ST)),SC),state="dimer",V1=nd$V,V2=nd$V)
dd[dd$T>4.5,"V2"]=dV[8:20,2]
dd[dd$C>4.5,"V2"]=dV[24:37,2]
#dd
dd=transform(dd,V=(V1+V2)/2)
dd=subset(dd,select=c(T,C,state,V))
dd

d=read.csv("barroso03tet.csv",header=T)
nd=rbind(data.frame(T=d[1:22,1],C=0,Vm=537,y=1/d[1:22,2]),
    data.frame(T=0,C=d[23:45,1],Vm=264,y=1/d[23:45,2]))
nd = transform(nd,V=Vm/(y+1))
dV=read.csv("barroso03Vtet.csv",header=T) # this is crap relative to log scale data
dV
dtet=data.frame(T=c(ST,rep(0,length(SC))),C=c(rep(0,length(ST)),SC),state="tetramer",V1=nd$V,V2=nd$V)
dtet[dtet$T>4.5,"V2"]=dV[11:23,2]
dtet[dtet$C>14,"V2"]=dV[30:39,2]
#dtet
dtet=transform(dtet,V=(V1+V2)/2)
dtet=subset(dtet,select=c(T,C,state,V))
dtet

d03=rbind(dd,dtet)
d03=data.frame(d03,year=2003)
d03

d=read.csv("barroso05dimer.csv",header=T)
#Vm=c(998,1261,463,732)
nd=rbind(data.frame(T=d[1:23,1],C=0,Vm=998,y=1/d[1:23,2]),
    data.frame(T=0,C=d[24:46,1],Vm=463,y=1/d[24:46,2]))
nd = transform(nd,V=Vm/(y+1))
nd
dV=read.csv("barroso05Vdimer.csv",header=T) # this is crap relative to log scale data
dV
ST=SC
dd=data.frame(T=c(ST,rep(0,length(SC))),C=c(rep(0,length(ST)),SC),state="dimer",V1=nd$V,V2=nd$V)
dd[dd$T>6,"V2"]=dV[12:24,2]
dd[dd$C>4.5,"V2"]=dV[33:46,2]
dd

dd=transform(dd,V=(V1+V2)/2)
dd=subset(dd,select=c(T,C,state,V))
dd

d=read.csv("barroso05tet.csv",header=T)
nd=rbind(data.frame(T=d[1:23,1],C=0,Vm=1261,y=1/d[1:23,2]),
    data.frame(T=0,C=d[24:46,1],Vm=732,y=1/d[24:46,2]))
nd = transform(nd,V=Vm/(y+1))
nd
dV=read.csv("barroso05Vtet.csv",header=T) # this is crap relative to log scale data
dV
dtet=data.frame(T=c(ST,rep(0,length(SC))),C=c(rep(0,length(ST)),SC),state="tetramer",V1=nd$V,V2=nd$V)
dtet[dtet$T>6,"V2"]=dV[12:24,2]
dtet[dtet$C>6,"V2"]=dV[34:46,2]
dtet
dtet=transform(dtet,V=(V1+V2)/2)
dtet=subset(dtet,select=c(T,C,state,V))
dtet

d05=rbind(dd,dtet)
d05=data.frame(d05,year=2005)
d05
write.csv(rbind(d03,d05),file="barroso.csv",row.names=FALSE)

