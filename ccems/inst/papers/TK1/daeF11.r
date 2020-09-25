if (.Platform$OS.type=="windows") home="/users/radivot" else 
  home="/home/radivot"
setwd(home)  # to place where models directory will go
library(ccems) 
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)  
g <-mkg(topology)
bK = c(0.85, 0.69, 0.65, 0.51)
kis = c(3.3, 3.9, 4.1, 4.1)
ST = 0.05
ET = 0.05 
names(bK)<-c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
names(kis)<-c("kE1S1","kE1S2","kE1S3","kE1S4")
p=c(bK,kis)
g$parmsTCC["ET"]=ET;
g$parmsTCC["ST"]=ST
g$parmsTCC["Kj_E1S1"]=p["E1S0_S"]/4  # convert to complete K
g$parmsTCC["Kj_E1S2"]=p["E1S0_S"]*p["E1S1_S"]/6
g$parmsTCC["Kj_E1S3"]=p["E1S0_S"]*p["E1S1_S"]*p["E1S2_S"]/4
g$parmsTCC["Kj_E1S4"]=p["E1S0_S"]*p["E1S1_S"]*p["E1S2_S"]*p["E1S3_S"]
dyn.load(strn<-paste("models/",g$id,.Platform$dynlib.ext,sep=""))
times <- seq(0,100,10)
out1 <- data.frame(lsoda(g$initialStateTCC,times,"myderivs", g$parmsTCC, 
        rtol=g$rtol,atol=g$atol, dllname=g$id))[11,2:3]
Eic=out1[1,"E"]  # the point of the code above was to get
Sic=out1[1,"S"]  # these three DAE initial conditions
ESic=ET-Eic

pars=c(p,g$parmsTCC)
frp=function(t,y,pars) {
  with (as.list(c(y, pars)), {
        E1S1 = S/Kj_E1S1
        E1S2 = S*S/Kj_E1S2
        E1S3 = S*S*S/Kj_E1S3
        E1S4 = S*S*S*S/Kj_E1S4
        denom=1+E1S1+E1S2+E1S3+E1S4
        num=(1/4)*(1*kE1S1*E1S1+2*kE1S2*E1S2+3*kE1S3*E1S3+4*kE1S4*E1S4)
        k=num/denom
        return(list(dy = c(dS=-4*k*ET, dP=4*k*ET)))
      })
}
dae <- function (t, y, yprime, pars) {
  with (as.list(c(y, yprime, pars)), {
        E1S1 = S/Kj_E1S1
        E1S2 = S*S/Kj_E1S2
        E1S3 = S*S*S/Kj_E1S3
        E1S4 = S*S*S*S/Kj_E1S4
        eq1 = ST - S - E*(E1S1+2*E1S2+3*E1S3+4*E1S4)      
        eq2 = ET - E - E*(E1S1+E1S2+E1S3+E1S4)
        res3= dST + E*(kE1S1*E1S1+2*kE1S2*E1S2+3*kE1S3*E1S3+4*kE1S4*E1S4)
        res4= dP + dST
        return(list(c(eq1, eq2, res3, res4)))  
      })
}
library(deSolve)
times <- seq(0, 10, by = .2)
y     <- c(S = ST, P = 0)
ODE <- as.data.frame(daspk(y = y, times = times, func = frp,
        parms = pars, atol = 1e-10, rtol = 1e-10))

y     <- c(ST = ST,S=Sic,E=Eic, P = 0)
dy  <- c(dST = 0, dS=0,dE=0, dP = 0) 
DAE <- as.data.frame(daspk(y = y, dy = dy, times = times,
        res = dae, parms = pars, atol = 1e-10, rtol = 1e-10))
par(mfrow = c(1,2))
for (i in c("S","P")) {
  plot(ODE$time,ODE[,i],xlab = "seconds", ylab = "micromolar",main = i,type = "l")
  lines(DAE$time,DAE[,i],col = "black",lty=2)
}
legend("bottomright",bty="n",lty = c(1,2), legend = c("ODE","DAE"))
