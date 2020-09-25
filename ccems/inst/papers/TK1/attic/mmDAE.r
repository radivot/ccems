rm(list=ls(all=TRUE))  # first clean up
if (.Platform$OS.type=="windows") 
             home="/users/radivot" else home="/home/radivot"
setwd(home)  # models directory (where C code lives) is off of home
library(ccems) 
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            m=c("E1S1")   
        ) # m for monomer 
    )
)  
g <-mkg(topology)
kmax=1;Km=1;ET=1;ST=2
g$parmsTCC["ET"]=ET;g$parmsTCC["ST"]=ST
g$parmsTCC["Kj_E1S1"]=Km
dyn.load(strn<-paste("models/",g$id,.Platform$dynlib.ext,sep=""))
times <- seq(0,100,10)
out1 <- data.frame(lsoda(g$initialStateTCC,times,"myderivs", g$parmsTCC, 
                   rtol=g$rtol,atol=g$atol, dllname=g$id))[11,2:3]
Eic=out1[1,"E"]
Sic=out1[1,"S"]
ESic=ET-Eic

vic=kmax*ESic
dSTdtic=-vic
del=0.1
g$parmsTCC["ST"]=g$parmsTCC["ST"]+del*dSTdtic
out2 <- data.frame(lsoda(g$initialStateTCC,times,"myderivs", g$parmsTCC, 
                  rtol=g$rtol,atol=g$atol, dllname=g$id))[11,2:3]
dFdtics=(out2-out1)/del
dSdtic=dFdtics[1,"S"]
dEdtic=dFdtics[1,"E"]

dyn.unload(strn)

pars <- c(kmax=kmax,Km=Km,ET=ET)

mmODE <- function (t, y, pars) {
  with (as.list(c(y, pars)), {
        k=kmax*S/Km/(1+S/Km)
        return(list(dy = c(dS=-k*ET, dP=k*ET)))
      })
}

mmDAE <- function (t, y, yprime, pars) {
  with (as.list(c(y, yprime, pars)), {
        res1 = dST - dS - dE*S/Km - E*dS/Km      
        res2 =      -dE - dE*S/Km - E*dS/Km
        ES=E*S/Km
        res3= dST + kmax*ES
        res4= dP + dST
        return(list(c(res1, res2, res3, res4),c(ES=ES,QST=S+ES)))  
      })
}

mmDAE2 <- function (t, y, yprime, pars) {
  with (as.list(c(y, yprime, pars)), {
        eq1 = ST - S - E*S/Km      
        eq2 = ET - E - E*S/Km
        ES=E*S/Km
        res3= dST + kmax*ES
        res4= dP + dST
        return(list(c(eq1, eq2, res3, res4),c(ES=ES,QST=S+ES)))  
      })
}



library(deSolve)
times <- seq(0, 10, by = del)
y     <- c(S = ST, P = 0)
ODE <- as.data.frame(daspk(y = y, times = times, func = mmODE,
        parms = pars, atol = 1e-10, rtol = 1e-10))

y     <- c(ST = ST,S=Sic,E=Eic, P = 0)
dy  <- c(dST = dSTdtic, dS=dSdtic,dE=dEdtic, dP = -dSTdtic) 
DAE <- as.data.frame(daspk(y = y, dy = dy, times = times,
        res = mmDAE, parms = pars, atol = 1e-10, rtol = 1e-10))
dy  <- c(dST = 0, dS=0,dE=0, dP = 0) 
DAE2 <- as.data.frame(daspk(y = y, dy = dy, times = times,
        res = mmDAE2, parms = pars, atol = 1e-10, rtol = 1e-10))
par(mfrow = c(1,2))
for (i in c("S","P")) {
  plot(ODE$time,ODE[,i],xlab = "time", ylab = "conc",main = i,type = "l")
#  if (i=="P") lines(ODE$time,((ODE[2,i]-ODE[1,i])/del)*ODE$time)
  lines(DAE2$time,DAE2[,i],col = "blue")
  lines(DAE$time,DAE[,i],col = "black",lty=2)
}
legend("bottomright",lty = c(1,2), legend = c("ODE","DAE"))
