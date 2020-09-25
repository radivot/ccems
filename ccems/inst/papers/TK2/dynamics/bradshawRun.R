library(odesolve)
#
#(* Mathematica 5 model of mitochondrial deoxynucleotide metabolism and DNA replication                          *)
#(* Author: Patrick Bradshaw   E-mail: bradshaw@vbi.vt.edu    10/14/2004                                         *)
#(* The model simulates transport of deoxynucleosides, dNDPs, and dNTPs into and out of a mitochondrion,          *)
#(*    the phosphorylation and dephosphorylation of nucleosides and nucleotides and synthesis of mtDNA from dNTPs *)
#(*****************************************************************************************************************)
#(*  Choose one of the 4 cell types you which to simulate. *)
#(*  cancer cells          =1                              *)
#(*  radidly dividing cells=2                               *)
#(*  slowly dividing cells =3                              *)
#(*  postmitotic cells     =4                               *)
celltype=3; # this switch is used within this constants file
source("C:\\Users\\radivot\\case\\active\\models\\dNTPsupply\\bradshaw\\bradshawModel.R");  # (* Read in constants file *)
y0=c(dT=dT0,dTMP=dTMP0,dTDP=dTDP0,dTTP=dTTP0, dC=dC0,dCMP=dCMP0,dCDP=dCDP0,dCTP=dCTP0,
    dG=dG0,dGMP=dGMP0,dGDP=dGDP0,dGTP=dGTP0, dA=dA0,dAMP=dAMP0,dADP=dADP0,dATP=dATP0, HDNA=HDNA0,LDNA=LDNA0)
y0
out1=lsoda(y=y0,times=seq(-50,0,1),fmito, parms=c(Hon=0,Lon=0), rtol=1e-4, atol= rep(1e-4,18)) 
yss=out1[nrow(out1),2:19]
yss

out2=lsoda(y=yss,times=seq(0,70,1),fmito, parms=c(Hon=1,Lon=0), rtol=1e-4, atol= rep(1e-4,18))
y0=out2[nrow(out2),2:19]
out3=lsoda(y=y0,times=seq(70,150,1),fmito, parms=c(Hon=0,Lon=0), rtol=1e-4, atol= rep(1e-4,18))
outs=data.frame(rbind(out1,out2,out3))

attach(outs)
par(mfrow=c(2,1))
plot(time,PolrateH,type="l")
plot(time,dTTP,type="l")

par(mfrow=c(1,1))
detach(outs)

